`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 24.06.2023 17:44:17
// Design Name: 
// Module Name: top_level_tb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module top_level_tb;

localparam NFFT = 256;
localparam NFFT_LOG2 =8; 
localparam BLOCO = 64;
localparam CHANNEL_SIZE = 64;
localparam NUM_CHANNELS = 6;

localparam SAMPLE_NUM = 283552;
localparam SAMPLE_NUM_LOG2 = 19;

localparam AWAIT_INIT_CONFIG = 0,
           RECEIVE_AUDIO = 1,
           FINISHED_AUDIO = 2,
           DONE = 3,
           SHIFT_TEMP_OUT = 4,
           ADD_TEMP_OUT = 5;

genvar i, j;
integer index, channel;

reg [2:0] tb_state;

reg clk;
reg rst;
reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] ts_audio;
reg sample_vad;

reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] audio_in [0:SAMPLE_NUM-1]; //era_apenas CHANNEL_SIZE
reg [SAMPLE_NUM_LOG2-1:0] sample_num, sample_counter;
reg [NFFT_LOG2-1:0] sample_addr;

reg sample_vads_in [0:SAMPLE_NUM-1];
reg ts_audio_valid = 0;
reg ts_audio_last = 0;
wire change_config;
integer addr = 0;
integer addr_out = 0;
reg fft_done = 0;

reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] mem [0:511];
reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] mem_out [0:511];

wire [CHANNEL_SIZE-1:0] input_data [0: NUM_CHANNELS-1];
wire [CHANNEL_SIZE/2-1:0] input_data_re [0: NUM_CHANNELS-1];
wire [CHANNEL_SIZE/2-1:0] input_data_im [0: NUM_CHANNELS-1];

wire [CHANNEL_SIZE-1:0] output_data [0: NUM_CHANNELS-1];
wire [CHANNEL_SIZE/2-1:0] output_data_re [0: NUM_CHANNELS-1];
wire [CHANNEL_SIZE/2-1:0] output_data_im [0: NUM_CHANNELS-1];
wire out_data_valid;
wire out_data_last;
wire [16*NFFT-1:0] out_data_info;

wire event_frame_started;
wire event_tlast_missing;
wire event_tlast_unexpected;
wire event_status_channel_halt;
wire event_data_in_channel_halt;
wire event_data_out_channel_halt;

wire config_complete;
wire fft_ready;
wire ready_for_sample;

wire [2*CHANNEL_SIZE-1:0] frame_out_packed;
reg [2*CHANNEL_SIZE-1:0] frame_out_temp [0:NFFT-1];
reg [2*CHANNEL_SIZE-1:0] frame_out_temp_to_sum [0:NFFT-1];
wire signed [CHANNEL_SIZE/2-1:0] frame_out_re [0:1];
wire signed [CHANNEL_SIZE/2-1:0] frame_out_im [0:1];
reg signed [CHANNEL_SIZE/2-1:0] frame_sum_out_re [0:NFFT-1][0:1];
reg signed [CHANNEL_SIZE/2-1:0] frame_sum_out_im [0:NFFT-1][0:1];
reg [CHANNEL_SIZE-1:0] final_audio_to_write [0:NFFT-1];
reg [7:0] temp_addr;

initial begin
    $readmemh("audio_input_fp_tb.txt", audio_in); 
    $readmemb("vad_in_tb.txt", sample_vads_in);
end

generate
    //for (i=0; i<NFFT; i=i+1) begin 
        for (j=0; j<2; j=j+1) begin
            assign frame_out_re[j] = frame_out_packed[j*CHANNEL_SIZE +: CHANNEL_SIZE/2];
            assign frame_out_im[j] = frame_out_packed[j*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
        end
    //end
endgenerate

top_level dut(
    .clk(clk),
    .rst(rst),
    .ts_audio_i(ts_audio),
    .ts_valid_i(ts_audio_valid),
    //.ts_last_i(ts_audio_last),
    .sample_vad(sample_vad),
    .event_fft_overflow_o(event_frame_started),
    .event_ifft_overflow_o(event_tlast_missing),
    .config_complete_o(config_complete),
    .fft_ready_o(fft_ready),
    .ready_for_sample_o(ready_for_sample),
    .ts_output_o(frame_out_packed),
    .ts_out_valid_o(out_data_valid),
    .ts_out_info(), //out_data_info
    //.ts_out_last_o(out_data_last)
    .ts_out_counter_o()
);


initial begin
    clk = 1'b0;
    forever #1 clk = ~clk;
end

initial begin
    rst = 1'b1;
    #6
    rst = 1'b0;
end

always @(posedge(clk)) begin
    if (rst) begin 
        ts_audio <= 0;
        fft_done <= 0;
        tb_state <= AWAIT_INIT_CONFIG;
        temp_addr <= 0;
        sample_addr <= 0;
        ts_audio_valid <= 0;
        sample_num <= 0;
        for (index=0; index<NFFT; index=index+1) begin
            for (channel=0; channel<2; channel=channel+1) begin
                frame_sum_out_re[index][channel] <= 0;
                frame_sum_out_im[index][channel] <= 0;
            end
        end
    end else begin
        case (tb_state)
            AWAIT_INIT_CONFIG: begin
                ts_audio_valid <= 0;
                sample_num <= 0;
                if (config_complete) begin
                    tb_state <= RECEIVE_AUDIO;
                    sample_counter <= 0;
                end
            end
            RECEIVE_AUDIO: begin
                if (ready_for_sample && fft_ready && sample_num < SAMPLE_NUM && sample_counter < BLOCO) begin 
                    ts_audio <= audio_in[sample_num];
                    sample_counter <= sample_counter + 1;
                    sample_vad <= sample_vads_in[sample_num];
                    ts_audio_valid <= 1;
                    sample_num <= sample_num + 1;
                end else begin
                    ts_audio_valid <= 0;
                end
                if (out_data_valid) begin 
                    frame_out_temp[temp_addr] <= frame_out_packed;
                    temp_addr <= temp_addr + 1;
                    if (temp_addr == NFFT-1) begin
                        tb_state <= SHIFT_TEMP_OUT;
                        temp_addr <= 0;
                    end
                end
            end
            FINISHED_AUDIO: begin
                
                
            
                //final_audio_to_write[sample_addr] <= {frame_sum_out_re[sample_addr][0], frame_sum_out_re[sample_addr][1]};
                final_audio_to_write[sample_addr] <= {frame_out_temp_to_sum[sample_addr][0 +: CHANNEL_SIZE/2], frame_out_temp_to_sum[sample_addr][CHANNEL_SIZE +: CHANNEL_SIZE/2]};
                sample_addr <= sample_addr + 1;
                if (sample_addr == NFFT-1) begin
                    tb_state <= DONE;
                    fft_done <= 1;
                end
            end
            SHIFT_TEMP_OUT: begin
                if (temp_addr+BLOCO < NFFT) begin
                    frame_out_temp_to_sum[temp_addr] <= frame_out_temp_to_sum[temp_addr+BLOCO];
                    temp_addr <= temp_addr + 1;
                end else begin
                    if (temp_addr < NFFT-1) begin 
                        frame_out_temp_to_sum[temp_addr] <= 0;
                        temp_addr <= temp_addr + 1;
                    end else begin
                        frame_out_temp_to_sum[temp_addr] <= 0;
                        tb_state <= ADD_TEMP_OUT;
                        temp_addr <= 0;
                    end
                end
            end
            ADD_TEMP_OUT: begin
                if (temp_addr < NFFT) begin
                    frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE/2-1:0] <= frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE/2-1:0] + frame_out_temp[temp_addr][CHANNEL_SIZE/2-1:0];
                    frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE-1:CHANNEL_SIZE/2] + frame_out_temp[temp_addr][CHANNEL_SIZE-1:CHANNEL_SIZE/2];
                    
                    frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE/2-1:CHANNEL_SIZE] <= frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE/2-1:CHANNEL_SIZE] + frame_out_temp[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE/2-1:CHANNEL_SIZE];
                    frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE-1:CHANNEL_SIZE/2+CHANNEL_SIZE] <= frame_out_temp_to_sum[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE-1:CHANNEL_SIZE/2+CHANNEL_SIZE] + frame_out_temp[temp_addr][CHANNEL_SIZE+CHANNEL_SIZE-1:CHANNEL_SIZE/2+CHANNEL_SIZE];
                    temp_addr <= temp_addr + 1;
                end 
                if (temp_addr == NFFT-1) begin
                    if (sample_num == SAMPLE_NUM) begin
                        tb_state <= FINISHED_AUDIO;
                    end else begin
                        sample_counter <= 0;
                        sample_addr <= 0;
                        tb_state <= RECEIVE_AUDIO;
                    end
                end
            end
            DONE: begin
            end
            default:
                tb_state <= AWAIT_INIT_CONFIG;
        endcase
        
    end
end

always @(fft_done) begin
    if (fft_done == 1) begin 
        $writememh("audio_output_fp_tb.txt", final_audio_to_write); 
        $finish;
    end
end

endmodule


