`timescale 1ns/1ps

module top_level #(
  //parameter SAMPLE_RATE = 441000,
  parameter NFFT = 256,
  parameter INTEGER_BITS = 12,
  parameter DECIMAL_BITS = 20,
  parameter CHANNEL_SIZE = 64,
  parameter NUM_CHANNELS = 6,
  parameter PRIOR_EST = 1001,
  parameter INPUT_WIDTH = CHANNEL_SIZE*NUM_CHANNELS,
  parameter Fs = 16000,
  parameter CONFIG_SIZE = 88,
  parameter TUSER_SIZE = 16 //log2(NFFT)+NUM_CHANNELS+PAD TO POWER OF 2 
)
(
  input wire clk,
  input wire rst,
  input wire [INPUT_WIDTH-1:0] ts_audio_i,
  input wire ts_valid_i,
  //input wire sample_window_i,
  //input wire ts_last_i,
  input wire sample_vad,
  output wire event_fft_overflow_o,
  output wire event_ifft_overflow_o,
  output wire config_complete_o,
  output wire fft_ready_o,
  output wire ready_for_sample_o,
  output wire [2*CHANNEL_SIZE-1:0] ts_output_o,
  output wire ts_out_valid_o,
  output wire [15:0] ts_out_info,
  output wire [31:0] ts_out_counter_o
);

localparam OUTPUT_WIDTH = 2*CHANNEL_SIZE;
localparam MULT_AMOUNT = 32;

localparam FWD = 1'b1;
localparam INV     = 1'b0;
localparam BLOCO = NFFT / 4;

//FSM states
localparam INIT = 0,
           ASSEMBLE_FRAME = 1,
           MULT_WEIGHT = 2,
           PROC_FRAME = 3,
           AWAIT_FREQ_MODULE = 4,
           IFFT = 5;

localparam MULT_TIME = 0,
           MULT_FREQ = 1;

genvar i, j;
integer index, channel;

reg [2:0] state;

reg signed [CHANNEL_SIZE/2-1:0] weight [0:NFFT-1];

wire core_ready;
wire freq_module_ready;
reg ready_for_sample_reg;
reg first_frame;
reg [NFFT-1:0] stored_vads;

reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] keep_bloco [0:BLOCO-1];

wire [CHANNEL_SIZE-1:0] input_channels [0:NUM_CHANNELS-1];

reg [INPUT_WIDTH-1:0] frame_buffer_1 [0: NFFT-1];
reg [INPUT_WIDTH-1:0] frame_buffer_2 [0: NFFT-1];

wire [2*NFFT*CHANNEL_SIZE-1:0] to_ifft_frame_packed;
reg [CHANNEL_SIZE-1:0] to_ifft_frame [0:NFFT-1][0:1];

reg [INPUT_WIDTH-1:0] fft_in_data_reg;
reg fft_in_valid_reg;
reg fft_in_last_reg, in_last_aux;

wire [INPUT_WIDTH-1:0] fft_in_data;
wire fft_in_valid;
wire fft_in_last;

wire [INPUT_WIDTH-1:0] freq_tdata;
wire freq_tvalid;
wire [TUSER_SIZE-1:0] freq_tuser;
wire freq_tlast;

reg [OUTPUT_WIDTH-1:0] to_ifft_tdata_reg;
reg to_ifft_tvalid_reg;
reg [TUSER_SIZE-1:0] freq_tuser_reg;
reg to_ifft_tlast_reg;

wire [OUTPUT_WIDTH-1:0] to_ifft_tdata;
wire to_ifft_tvalid;
wire to_ifft_tlast;

reg config_complete;
reg fft_config_complete;
reg ifft_config_complete;

reg [CHANNEL_SIZE-1:0] output_channels [0:NUM_CHANNELS-1];

reg [87:0] fft_config_tdata_reg;
reg fft_config_tvalid_reg;

wire [87:0] fft_config_tdata;
wire fft_config_tvalid;
wire fft_config_tready;

reg [23:0] ifft_config_tdata_reg;
reg ifft_config_tvalid_reg;

wire [23:0] ifft_config_tdata;
wire ifft_config_tvalid;
wire ifft_config_tready;

wire ifft_ready;

wire [7:0] fft_status_tdata;
wire fft_status_tvalid;
wire fft_overflow;

wire [7:0] ifft_status_tdata;
wire ifft_status_tvalid;
wire ifft_overflow;

reg [7:0] sample_counter;
reg [7:0] sample_counter2;
reg buffer_choice;

reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] assembled_frame [0:NFFT-1];
reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] assembled_frame_after_weight [0:NFFT-1];
reg [7:0] assembled_addr;
reg [7:0] to_freq_addr;
//reg [NFFT-1:0] to_ifft_addr;
//reg [NFFT-1:0] after_ifft_addr;
reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] to_freq_frame [0:NFFT-1];
//wire [CHANNEL_SIZE*NUM_CHANNELS*NFFT-1:0] to_freq_frame_packed;
//wire [CHANNEL_SIZE*NUM_CHANNELS*NFFT-1:0] to_ifft_frame_packed;
//reg [CHANNEL_SIZE-1:0] to_ifft_frame [0:NUM_CHANNELS-1][0:NFFT-1];
wire to_ifft_frame_valid;
reg to_freq_frame_valid;
reg signed [CHANNEL_SIZE/2-1:0] assembled_frame_to_mult_re [0:NFFT-1][0:NUM_CHANNELS-1];
reg signed [CHANNEL_SIZE/2-1:0] assembled_frame_to_mult_im [0:NFFT-1][0:NUM_CHANNELS-1];

//reg signed [CHANNEL_SIZE/2-1:0] after_mult_frame_to_fft_re [0:NUM_CHANNELS-1][0:NFFT-1];
//reg signed [CHANNEL_SIZE/2-1:0] after_mult_frame_to_fft_im [0:NUM_CHANNELS-1][0:NFFT-1];
reg signed [CHANNEL_SIZE-1:0] after_mult_frame_to_fft [0:NFFT-1][0:NUM_CHANNELS-1];
reg after_weight_frame_valid;

reg [CHANNEL_SIZE-1:0] in_a [0:NFFT-1][0:NUM_CHANNELS-1];
reg [CHANNEL_SIZE-1:0] in_b [0:NFFT-1][0:NUM_CHANNELS-1];
wire [CHANNEL_SIZE-1:0] result [0:NFFT-1][0:NUM_CHANNELS-1];
//reg [CHANNEL_SIZE-1:0] result [0:NFFT-1][0:NUM_CHANNELS-1];
reg valid [0:NFFT-1][0:NUM_CHANNELS-1];
wire valid_w [0:NFFT-1][0:NUM_CHANNELS-1];
wire done [0:NFFT-1][0:NUM_CHANNELS-1];
wire ready [0:NFFT-1][0:NUM_CHANNELS-1];
//reg ready [0:NFFT-1][0:NUM_CHANNELS-1];
wire ovflw_w [0:NFFT-1][0:NUM_CHANNELS-1];

reg [1:0] freq_phase_pattern_reg;
wire freq_total_vad;
reg [0:NFFT-1] freq_sample_vads;

wire [2*CHANNEL_SIZE-1:0] ifft_core_out_tdata;
wire [15:0] ifft_core_out_tuser;
wire ifft_core_out_tvalid;
wire ifft_core_out_tlast;

reg [2*CHANNEL_SIZE-1:0] after_ifft_frame [0:NFFT-1];
reg after_ifft_frame_valid;

reg [CHANNEL_SIZE/2-1:0] in_a_re [0:NFFT-1][0:NUM_CHANNELS-1];
reg [CHANNEL_SIZE/2-1:0] in_b_re [0:NFFT-1][0:NUM_CHANNELS-1];
wire [CHANNEL_SIZE/2-1:0] result_re [0:NFFT-1][0:NUM_CHANNELS-1];
//reg [CHANNEL_SIZE/2-1:0] result_re [0:NFFT-1][0:NUM_CHANNELS-1];
reg valid_re [0:NFFT-1][0:NUM_CHANNELS-1];
wire done_re [0:NFFT-1][0:NUM_CHANNELS-1];
wire ready_re [0:NFFT-1][0:NUM_CHANNELS-1];
wire ovflw_re [0:NFFT-1][0:NUM_CHANNELS-1];
//reg done_re [0:NFFT-1][0:NUM_CHANNELS-1];
//reg ready_re [0:NFFT-1][0:NUM_CHANNELS-1];
//reg ovflw_re [0:NFFT-1][0:NUM_CHANNELS-1];

reg [CHANNEL_SIZE/2-1:0] final_frame_out_re [0:NFFT-1][0:1];
reg [2*CHANNEL_SIZE-1:0] final_word_out;
wire [2*CHANNEL_SIZE*NFFT-1:0] final_frame_out_packed;

reg final_frame_valid_reg;
reg [16*NFFT-1:0] ts_out_info_reg;
reg [31:0] frame_counter_reg;

//MULT SIGNALS
reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_in_a;
reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_in_b;
reg [MULT_AMOUNT-1:0] mult_in_valid;
wire [MULT_AMOUNT-1:0] mult_in_real_mult;
wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_out_result;
wire [MULT_AMOUNT-1:0] mult_out_done;
wire [MULT_AMOUNT-1:0] mult_out_ready;
wire [MULT_AMOUNT-1:0] mult_out_ovflw;

reg [MULT_AMOUNT*CHANNEL_SIZE-1:0]  mult_a_time, mult_b_time;
wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_a_freq, mult_b_freq;
reg [MULT_AMOUNT-1:0]  mult_time_valid, freq_mults_done;
wire [MULT_AMOUNT-1:0] mult_freq_valid;
reg to_freq_vad, bloco_vad, freq_word_valid;
reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] assembled_frame_dout, assembled_frame_din, freq_frame_dout;
reg [CHANNEL_SIZE/2-1:0] weight_dout;
reg [7:0] addr, ifft_addr_in, ifft_addr_out;
reg [2*CHANNEL_SIZE-1:0] ifft_frame_dout;

wire ready_ifft;
reg after_ifft_mult_en;
reg [1:0] ifft_mults_done;
wire to_ifft_tvalid_corr, to_ifft_tlast_corr;
wire [2*CHANNEL_SIZE-1:0] to_ifft_tdata_corr;
reg check_empty_re;// check_empty_im;
reg mult_mux, back_to_assemble;
wire fft_tlast_missing, fft_tlast_unexp;

assign ts_output_o = final_word_out;

assign mult_in_real_mult = 0;

initial begin 
    $readmemh("weight_fp_hex.txt", weight);
    $readmemh("zero_init_mem.txt", assembled_frame);
end

assign fft_in_valid = fft_in_valid_reg;
assign fft_in_last = fft_in_last_reg;

assign freq_total_vad = | freq_sample_vads;
assign ts_out_info = ts_out_info_reg;

assign fft_in_data = fft_in_data_reg;
//CONFIG DATA => | PAD | SCALE_SCH | FWD/INV | PAD | CP_LEN | PAD | NFFT |

assign config_complete_o = config_complete;

assign fft_config_tdata = fft_config_tdata_reg;
assign fft_config_tvalid = fft_config_tvalid_reg;



assign fft_ready_o = core_ready && config_complete;
assign ready_for_sample_o = ready_for_sample_reg;

multiplicators #(
    .DECIMAL_BITS(DECIMAL_BITS),
    .CHANNEL_SIZE(CHANNEL_SIZE),
    .MULT_AMOUNT(MULT_AMOUNT)
) mult_inst (
    .clk(clk),
    .rst(rst),
    .input_a(mult_in_a),
    .input_b(mult_in_b),
    .valid(mult_in_valid),
    .real_mult(mult_in_real_mult),
    .result(mult_out_result),
    .ready(mult_out_ready),
    .done(mult_out_done),
    .ovflw(mult_out_ovflw)
);

xfft_0 fft (
  .aclk(clk),                                                // input wire aclk
  .aresetn(~rst),                                          // input wire aresetn
  .s_axis_config_tdata(fft_config_tdata),                  // input wire [87 : 0] s_axis_config_tdata
  .s_axis_config_tvalid(fft_config_tvalid),                // input wire s_axis_config_tvalid
  .s_axis_config_tready(fft_config_tready),                // output wire s_axis_config_tready
  .s_axis_data_tdata(fft_in_data),                      // input wire [383 : 0] s_axis_data_tdata
  .s_axis_data_tvalid(fft_in_valid),                    // input wire s_axis_data_tvalid
  .s_axis_data_tready(core_ready),                    // output wire s_axis_data_tready
  .s_axis_data_tlast(fft_in_last),                      // input wire s_axis_data_tlast
  .m_axis_data_tdata(freq_tdata),                      // output wire [383 : 0] m_axis_data_tdata
  .m_axis_data_tuser(freq_tuser),                      // output wire [15 : 0] m_axis_data_tuser
  .m_axis_data_tvalid(freq_tvalid),                    // output wire m_axis_data_tvalid
  .m_axis_data_tready(1'b1),                    // input wire m_axis_data_tready
  .m_axis_data_tlast(freq_tlast),                      // output wire m_axis_data_tlast
  .m_axis_status_tdata(fft_status_tdata),                  // output wire [7 : 0] m_axis_status_tdata
  .m_axis_status_tvalid(fft_status_tvalid),                // output wire m_axis_status_tvalid
  .m_axis_status_tready(1'b1),                // input wire m_axis_status_tready
  .event_frame_started(),                  // output wire event_frame_started
  .event_tlast_unexpected(fft_tlast_unexp),            // output wire event_tlast_unexpected
  .event_tlast_missing(fft_tlast_missing),                  // output wire event_tlast_missing
  .event_fft_overflow(fft_overflow),                    // output wire event_fft_overflow
  .event_status_channel_halt(),      // output wire event_status_channel_halt
  .event_data_in_channel_halt(),    // output wire event_data_in_channel_halt
  .event_data_out_channel_halt()  // output wire event_data_out_channel_halt
);

assign event_fft_overflow_o = fft_status_tvalid ? fft_status_tdata : 8'h00;

assign to_ifft_tdata = to_ifft_tdata_reg;
assign to_ifft_tvalid = to_ifft_tvalid_reg;
assign to_ifft_tlast = to_ifft_tlast_reg;
assign ts_out_counter_o = frame_counter_reg;

assign ts_out_valid_o = final_frame_valid_reg;

frequency_proc_alt #(
    .NFFT(NFFT),
    .INTEGER_BITS(INTEGER_BITS),
    .DECIMAL_BITS(DECIMAL_BITS),
    .CHANNEL_SIZE(CHANNEL_SIZE),
    .NUM_CHANNELS(NUM_CHANNELS),
    .PRIOR_EST(PRIOR_EST)
) freq_inst(
    .clk(clk),
    .rst(rst),
    .vad_i(to_freq_vad),
    .freq_word_i(freq_frame_dout),  //mudar para uma palavra por vez
    .freq_wrd_valid_i(freq_word_valid),
    .frame_valid_i(to_freq_frame_valid),
    //.freq_tvalid_i(to_freq_frame_valid),
    .phase_pattern_i(freq_phase_pattern_reg),
    .mult_done(mult_out_done),
    .mult_ready(mult_out_ready),
    .mult_result(mult_out_result),
    .check_empty(check_empty_re),
    .to_mult_valid(mult_freq_valid),
    .to_mult_real_mult(),
    .to_mult_in_a(mult_a_freq),
    .to_mult_in_b(mult_b_freq),
    .ifft_core_ready_i(ifft_ready),
    .to_ifft_tdata_o(to_ifft_tdata_corr), //to_ifft_frame_packed
    .ready_ifft_o(ready_ifft),
    .to_ifft_tlast_o(to_ifft_tlast_corr),
    .ready_o(freq_module_ready),
    .to_ifft_tvalid_o(to_ifft_tvalid_corr) //to_ifft_frame_valid
);

xfft_1 ifft (
  .aclk(clk),                                                // input wire aclk
  .aresetn(~rst),                                          // input wire aresetn
  .s_axis_config_tdata(ifft_config_tdata),                  // input wire [23 : 0] s_axis_config_tdata
  .s_axis_config_tvalid(ifft_config_tvalid),                // input wire s_axis_config_tvalid
  .s_axis_config_tready(ifft_config_tready),                // output wire s_axis_config_tready
  .s_axis_data_tdata(to_ifft_tdata_corr),                      // input wire [127 : 0] s_axis_data_tdata
  .s_axis_data_tvalid(to_ifft_tvalid_corr),                    // input wire s_axis_data_tvalid
  .s_axis_data_tready(ifft_ready),                    // output wire s_axis_data_tready
  .s_axis_data_tlast(to_ifft_tlast_corr),                      // input wire s_axis_data_tlast
  .m_axis_data_tdata(ifft_core_out_tdata),                      // output wire [127 : 0] m_axis_data_tdata
  .m_axis_data_tuser(ifft_core_out_tuser),                      // output wire [15 : 0] m_axis_data_tuser
  .m_axis_data_tvalid(ifft_core_out_tvalid),                    // output wire m_axis_data_tvalid
  .m_axis_data_tready(1'b1),                    // input wire m_axis_data_tready
  .m_axis_data_tlast(ifft_core_out_tlast),                      // output wire m_axis_data_tlast
  .m_axis_status_tdata(ifft_status_tdata),                  // output wire [7 : 0] m_axis_status_tdata
  .m_axis_status_tvalid(ifft_status_tvalid),                // output wire m_axis_status_tvalid
  .m_axis_status_tready(1'b1),                // input wire m_axis_status_tready
  .event_frame_started(),                  // output wire event_frame_started
  .event_tlast_unexpected(),            // output wire event_tlast_unexpected
  .event_tlast_missing(),                  // output wire event_tlast_missing
  .event_fft_overflow(ifft_overflow),                    // output wire event_fft_overflow
  .event_status_channel_halt(),      // output wire event_status_channel_halt
  .event_data_in_channel_halt(),    // output wire event_data_in_channel_halt
  .event_data_out_channel_halt()  // output wire event_data_out_channel_halt
);
assign ifft_config_tdata = ifft_config_tdata_reg;
assign ifft_config_tvalid = ifft_config_tvalid_reg;
assign event_ifft_overflow_o = ifft_status_tvalid ? ifft_status_tdata : 8'h00;

always @(posedge(clk)) begin
    to_ifft_tlast_reg <= 0;
    fft_in_valid_reg <= 0;
    to_ifft_tvalid_reg <= 0;
    final_frame_valid_reg <= 0;
    fft_in_last_reg <= 0;
    to_freq_frame_valid <= 0;
    if (rst) begin
        check_empty_re <= 0;
        back_to_assemble <= 0;
        //check_empty_im <= 0;
        mult_a_time <= 0;
        mult_b_time <= 0;
        mult_mux <= MULT_TIME;
        stored_vads <= 0;
        freq_phase_pattern_reg <= 2'b11;
        freq_sample_vads <= 0;
        frame_counter_reg <= 0;
        freq_word_valid <= 0;
        //addr <= 0;
        state <= INIT;
        valid_re[index][channel] <= 0;
        config_complete <= 1'b0;
        ready_for_sample_reg <= 0;
        after_weight_frame_valid <= 0;
        to_freq_vad <= 0;
        bloco_vad <= 0;
        first_frame <= 1;
        assembled_addr <= 0;
        fft_config_complete <= 1'b0;
        ifft_config_complete <= 1'b0;
        fft_config_tdata_reg[5:0] = {6{FWD}}; //FWD/INV
        fft_config_tdata_reg[53:6] = {48{1'b0}}; //SCALING
        fft_config_tdata_reg[87:54] = {34{1'b0}}; //PAD
        fft_config_tvalid_reg <= 1'b1;
        
        ifft_config_tdata_reg <= {24{1'b0}};
        ifft_config_tvalid_reg <= 1'b1;
        
        to_ifft_tvalid_reg <= 1'b0;
        for (index = 0; index<NUM_CHANNELS; index=index+1) begin
            output_channels[index] <= {CHANNEL_SIZE{1'b0}};
        end
    end else begin
        mult_time_valid <= 0;
        case(state)
            INIT: begin
                //fft_config_tvalid_reg <= 1'b1;
                //ifft_config_tvalid_reg <= 1'b1;
                if (fft_config_tready && fft_config_tvalid_reg) begin 
                    fft_config_tvalid_reg <= 1'b0;
                    fft_config_complete <= 1'b1;
                end
                if (ifft_config_tready && ifft_config_tvalid_reg) begin 
                    ifft_config_tvalid_reg <= 1'b0;
                    ifft_config_complete <= 1'b1;
                end
                if (fft_config_complete && ifft_config_complete) begin 
                    state <= ASSEMBLE_FRAME; 
                    //stored_vads <= |stored_vads
                    //to_freq_vad <= |stored_vads[NFFT-BLOCO-1:0];
                    if (first_frame == 0) begin 
                        to_freq_vad <= 0;
                        assembled_addr <= 0;
                        ready_for_sample_reg <= 0;
                    end else begin 
                        assembled_addr <= NFFT-BLOCO;
                        ready_for_sample_reg <= 1;
                    end
                    
                    config_complete <= 1'b1;
                end
            end
            ASSEMBLE_FRAME:  begin
                if (first_frame == 0) begin 
                        //assembled_frame[assembled_addr] <= keep_bloco[assembled_addr];
                        
                        if (assembled_addr < NFFT-BLOCO) begin
                            assembled_frame[assembled_addr] <= assembled_frame[assembled_addr+BLOCO];
                            assembled_addr <= assembled_addr + 1; 
                            //first_frame <= 1;
                            
                            //to_freq_vad <= bloco_vad;
                        end else begin
                            ready_for_sample_reg <= 1;
                        end
//                        for (index=0; index<BLOCO; index=index+1) begin 
//                            for (channel=0; channel<NUM_CHANNELS; channel=channel+1) begin
//                                assembled_frame[index][channel] <=  keep_bloco[index][channel];
//                            end
//                        end
                end
                if (ts_valid_i && ready_for_sample_reg) begin
                    to_freq_vad <= to_freq_vad | sample_vad;
                    assembled_frame[assembled_addr] <= ts_audio_i;
                    assembled_addr <= assembled_addr + 1;
                    if (assembled_addr == NFFT-1) begin 
                        first_frame <= 0;
                        ready_for_sample_reg <= 0;
                        after_weight_frame_valid <= 0;
                        assembled_addr <= 0;
                        to_freq_addr <= 0;
                        to_freq_frame_valid <= 0;
                        assembled_frame_dout <= assembled_frame[0];
                        weight_dout <= weight[0];
                        state <= MULT_WEIGHT;
                        freq_mults_done <= 0;
                        if (freq_phase_pattern_reg == 2'b11) begin 
                            freq_phase_pattern_reg <= 0;
                        end else begin 
                            freq_phase_pattern_reg <= freq_phase_pattern_reg + 1;
                        end
                    end
                end
            end
            MULT_WEIGHT: begin
                for (index=0; index<NUM_CHANNELS; index=index+1) begin
                    if (mult_out_ready[index]) begin
                        mult_time_valid[index] <= 1;
                        mult_a_time[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= assembled_frame_dout[index*CHANNEL_SIZE +: CHANNEL_SIZE];
                        mult_b_time[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, weight_dout};
                    end
                    if (mult_out_done[index]) begin
                        freq_mults_done[index] <= 1;
                        assembled_frame_din[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= mult_out_result[index*CHANNEL_SIZE +: CHANNEL_SIZE];
                    end
                end
                if (&(freq_mults_done[NUM_CHANNELS-1:0])) begin
                    freq_mults_done <= 0;
                    assembled_frame_after_weight[assembled_addr] <= assembled_frame_din;
                    if (assembled_addr == NFFT-1) begin
                        assembled_addr <= 1;
                        after_weight_frame_valid <= 1;
                        check_empty_re <= 0;
                        state <= PROC_FRAME;
                        in_last_aux <= 0;
                        mult_mux <= MULT_FREQ;
                        assembled_frame_dout <= assembled_frame_after_weight[0];
                    end else begin
                        after_weight_frame_valid <= 0;
                        assembled_frame_dout <= assembled_frame[assembled_addr+1];
                        weight_dout <= weight[assembled_addr+1];
                        assembled_addr <= assembled_addr + 1;
                    end
                    
                end
            end
            PROC_FRAME: begin
                if (after_weight_frame_valid && assembled_addr > 0 && assembled_addr < NFFT) begin 
                    if (core_ready) begin 
                        assembled_frame_dout <= assembled_frame_after_weight[assembled_addr];
                        if (assembled_addr < NFFT-1) begin
                            assembled_addr <= assembled_addr + 1;
                        end
                        if (assembled_addr == NFFT-1) begin 
                            in_last_aux <= 1;
                            
                        end
                        if (in_last_aux) begin
                            fft_in_last_reg <= 1;
                            assembled_addr <= 0;
                        end else begin
                            fft_in_last_reg <= 0;
                        end
                        fft_in_data_reg <= assembled_frame_dout;
                        fft_in_valid_reg <= 1;
                    end
                end else begin 
                    fft_in_valid_reg <= 0;
                    fft_in_last_reg <= 0;
                end
                if (freq_tvalid && to_freq_addr < NFFT) begin 
                    check_empty_re <= check_empty_re | (|freq_tdata);
                    to_freq_frame[to_freq_addr] <= freq_tdata;
                    to_freq_addr <= to_freq_addr + 1;
                end
                if (to_freq_addr == NFFT-1 && freq_module_ready) begin 
                    to_freq_frame_valid <= 1;
                    //to_ifft_addr <= 0;
                    to_freq_addr <= 1;
                    //after_ifft_addr <= 0;
                    after_ifft_frame_valid <= 0;
                    
                    freq_frame_dout <= to_freq_frame[0];
                    if (freq_module_ready) begin
                        freq_word_valid <= 1;
                        state <= AWAIT_FREQ_MODULE;
                    end
                end
            end
            AWAIT_FREQ_MODULE: begin
                if (freq_module_ready) begin
                    if (to_freq_addr > 0 && to_freq_addr < NFFT) begin
                        freq_word_valid <= 1;
                        freq_frame_dout <= to_freq_frame[to_freq_addr];
                        to_freq_addr <= to_freq_addr + 1;
                    end else begin
                        freq_word_valid <= 0;
                    end
                end
                if (ready_ifft) begin
                    state <= IFFT;
                    mult_mux <= MULT_TIME;
                    after_ifft_mult_en <= 0;
                    ifft_addr_in <= 0;
                    ifft_mults_done <= 0;
                    ifft_addr_out <= 0;
                    //weight_dout <= weight[0];
                end
            end
            IFFT: begin
                mult_time_valid[1:0] <= 0;
                final_frame_valid_reg <= 0;
                if (ifft_core_out_tvalid && ifft_addr_in < NFFT) begin 
                    after_ifft_frame[ifft_addr_in] <= ifft_core_out_tdata;
                    if (ifft_addr_in < NFFT-1) begin
                        ifft_addr_in <= ifft_addr_in + 1;
                    end
                    
                    ts_out_info_reg[ifft_addr_in] <= ifft_core_out_tuser;
                end
                if (ifft_addr_in > 1 && after_ifft_mult_en == 0) begin
                        after_ifft_mult_en <= 1;
                        weight_dout <= weight[ifft_addr_out];
                        ifft_frame_dout <= after_ifft_frame[ifft_addr_out];
                        if (ifft_addr_out == NFFT-1) begin
                            back_to_assemble <= 1;
                        end
                        ifft_addr_out <= ifft_addr_out + 1;
                end
                for (index=0; index<2; index=index+1) begin
                    if (after_ifft_mult_en && mult_out_ready[index]) begin
                    mult_time_valid[index] <= 1;
                    mult_a_time[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= ifft_frame_dout[index*CHANNEL_SIZE +: CHANNEL_SIZE]; 
                    mult_b_time[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, weight_dout};
                    end
                    if (after_ifft_mult_en && mult_out_done[index]) begin
                        final_word_out[index*CHANNEL_SIZE +: CHANNEL_SIZE] <= mult_out_result[index*CHANNEL_SIZE +: CHANNEL_SIZE];
                        ifft_mults_done[index] <= 1;
                    end
                end
                if (&ifft_mults_done) begin
                    ifft_mults_done <= 0;
                    after_ifft_mult_en <= 0;
                    final_frame_valid_reg <= 1;
                    if (back_to_assemble) begin
                        back_to_assemble <= 0;
                        state <= ASSEMBLE_FRAME;
                    end
                end
            end
        endcase
    end
end

always @(mult_mux, mult_time_valid, mult_freq_valid) begin
    case(mult_mux)
        MULT_TIME: begin
            mult_in_valid = mult_time_valid;
        end
        MULT_FREQ: begin
            mult_in_valid = mult_freq_valid;
        end
    endcase
end
always @(mult_mux, mult_a_time, mult_a_freq) begin
    case(mult_mux)
        MULT_TIME: begin
            mult_in_a = mult_a_time;
        end
        MULT_FREQ: begin
            mult_in_a = mult_a_freq;
        end
    endcase
end
always @(mult_mux, mult_b_time, mult_b_freq) begin
    case(mult_mux)
        MULT_TIME: begin
            mult_in_b = mult_b_time;
        end
        MULT_FREQ: begin
            mult_in_b = mult_b_freq;
        end
    endcase
end


endmodule
