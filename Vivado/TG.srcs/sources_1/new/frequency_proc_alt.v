`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 23.08.2023 12:00:33
// Design Name: 
// Module Name: frequency_proc_alt
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


module frequency_proc_alt #(
    parameter NFFT = 256,
    parameter INTEGER_BITS = 12,
    parameter NUM_PATTERNS = 4,
    parameter DECIMAL_BITS = 20,
    parameter CHANNEL_SIZE = 64,
    parameter NUM_CHANNELS = 6,
    parameter MULT_AMOUNT = 32,
    parameter PRIOR_EST = 4000,
    parameter OUTPUT_WIDTH = 2*CHANNEL_SIZE,
    parameter INPUT_WIDTH = CHANNEL_SIZE*NUM_CHANNELS
)(
    input wire clk,
    input wire rst,
    input wire vad_i,
    input wire [INPUT_WIDTH-1:0] freq_word_i,
    input wire freq_wrd_valid_i,
    input wire frame_valid_i,
    input wire [1:0] phase_pattern_i,
    input wire [MULT_AMOUNT-1:0] mult_done,
    input wire [MULT_AMOUNT-1:0] mult_ready,
    input wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_result,
    input wire check_empty,
    output wire [MULT_AMOUNT-1:0] to_mult_valid,
    output wire [MULT_AMOUNT-1:0] to_mult_real_mult,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] to_mult_in_a,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] to_mult_in_b,
    output wire ready_o,
    input wire ifft_core_ready_i,
    output wire [OUTPUT_WIDTH-1:0] to_ifft_tdata_o,
    output wire ready_ifft_o,
    output wire to_ifft_tlast_o,
    output wire to_ifft_tvalid_o
    //output to_ifft_tlast_o
    );
    integer i;
    assign to_mult_real_mult = 0;
    localparam AWAIT_FRAME  = 0,
               PHASE_ADJUST = 1,
               AWAIT_CORR   = 2,
               AWAIT_MWF    = 3,
               AWAIT_CONJ_COEF_MULT   = 4,
               STORE_PHASED = 5,
               READY_FOR_IFFT = 6,
               TO_IFFT_CONJ_PHASE_MULT = 7,
               AWAIT_qLR_MULT = 8;
   
   localparam MULT_FREQ = 0,
              MULT_CORR = 1,
              MULT_MWF = 2;
              
   localparam CONC_MULTS = 5;
   localparam REDUCED_MULTS = 30;
   localparam MATRIX_WORD_SIZE = CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS;
   reg [3:0] state;
   reg [1:0] mult_sel;
   reg ready_reg;
   reg [31:0] cont_frame;
   
   reg [CHANNEL_SIZE-1:0] phase_0_mem [0:NFFT-1];  
   reg [CHANNEL_SIZE-1:0] phase_1_mem [0:NFFT-1];  
   reg [CHANNEL_SIZE-1:0] phase_2_mem [0:NFFT-1];  
   reg [CHANNEL_SIZE-1:0] phase_3_mem [0:NFFT-1]; 
   
   reg [NUM_CHANNELS*CHANNEL_SIZE/2-1:0] qL;
   reg [NUM_CHANNELS*CHANNEL_SIZE/2-1:0] qR;
   
   reg [INPUT_WIDTH-1:0] freq_frame [0:NFFT-1];
   reg [7:0] sample_addr, phase_addr;
   reg [2:0] phase_addr_aux, words_aux;
   reg [INPUT_WIDTH-1:0] frame_word_in, frame_word_out;
   reg [CONC_MULTS*INPUT_WIDTH-1:0] freq_word_to_phase, phase_to_mult_arr;
   reg [CONC_MULTS*NUM_CHANNELS-1:0] freq_to_phase_valid;
   reg frame_vad, word_valid, frame_prep_done_corr, frame_prep_done_mwf;
   
   wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_in_a_proc, mult_in_b_proc, mult_in_a_corr, mult_in_b_corr, mult_in_a_mwf, mult_in_b_mwf;
   reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] mult_in_a_proc_reg, mult_in_b_proc_reg, mult_in_a_selected, mult_in_b_selected;
   reg [CONC_MULTS*NUM_CHANNELS*CHANNEL_SIZE-1:0] stored_phase_res;
   reg[CONC_MULTS-1:0] finished_mults;
   wire [MULT_AMOUNT-1:0] mult_valid_proc, mult_out_ready, mult_out_done, mult_valid_corr, mult_valid_mwf;
   reg [MULT_AMOUNT-1:0] mult_valid_proc_reg, mult_valid_selected;
   reg [CHANNEL_SIZE-1:0] phase_to_mult;
   reg [4:0] words_read;
   
   reg matrix_corr_do_div;
   reg matrix_corr_est_mode, to_corr_ack;
   wire matrix_corr_done, from_corr_req;
   wire matrix_corr_ready;
   wire matrix_corr_div_done;
   wire [7:0] to_corr_address, to_mwf_address;
   
   wire [MATRIX_WORD_SIZE-1:0] corr_to_mwf_ry_word, corr_to_mwf_rx_word;
   reg mwf_matrix_w_updated, to_mwf_ack, from_mwf_req, m_corr_done;
   wire mwf_w_update_req, mwf_done, mwf_ready, mwf_word_updated;
   
   reg mwf_coef_req, phase_to_mult_reg;
   wire mwf_coef_ack;
   wire [INPUT_WIDTH*2-1:0] coef_sample;
   
   reg [2*CHANNEL_SIZE-1:0] ifft_data_after_coeff, to_ifft_tdata_reg, ifft_data_after_coeff_sum;
   reg [2*NUM_CHANNELS*CHANNEL_SIZE-1:0] ifft_data_after_coeff_to_sum;
   reg [2*NUM_CHANNELS-1:0] after_coeff_mult_done, qLR_mult_done, mult_remain;
   reg [CONC_MULTS*NUM_CHANNELS-1:0] mult_phase_remain;
   
   reg [2*CHANNEL_SIZE-1:0] to_ifft_mem [0:NFFT-1];
   reg ready_ifft_reg, to_ifft_tvalid_reg, to_ifft_tlast_reg;
   //reg [OUTPUT_WIDTH-1:0] to_ifft_tdata_reg;
   
   `ifdef SIMULATION
        localparam DESIRED_FRAME = 2;
        reg [31:0] frame_cont_debug;
        reg pre_ph_done, post_ph_done, post_qlr_done;
        reg [7:0] addr_pre_phase;
        reg [7:0] addr_post_phase;
        reg [7:0] addr_pre_qlr_phase;
        reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] pre_phase_debug [0:NFFT-1];
        reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] post_phase_debug [0:NFFT-1];
        reg [CHANNEL_SIZE*NUM_CHANNELS-1:0] post_QLR_debug [0:NFFT-1];
        
   `endif
   
   initial begin 
        $readmemh("phase_0.txt", phase_0_mem);
        $readmemh("phase_1.txt", phase_1_mem);
        $readmemh("phase_2.txt", phase_2_mem);
        $readmemh("phase_3.txt", phase_3_mem);
    end
   
   assign to_mult_in_a = mult_in_a_selected;
   assign to_mult_in_b = mult_in_b_selected;
   assign to_mult_valid = mult_valid_selected;
   
   assign mult_valid_proc = mult_valid_proc_reg;
   assign mult_in_a_proc = mult_in_a_proc_reg;
   assign mult_in_b_proc = mult_in_b_proc_reg;
   assign ready_ifft_o = ready_ifft_reg;
   assign to_ifft_tdata_o = to_ifft_tdata_reg;
   assign to_ifft_tvalid_o = to_ifft_tvalid_reg;
   assign to_ifft_tlast_o = to_ifft_tlast_reg;
   assign ready_o = ready_reg;
   
   always @(posedge(clk)) begin
    if (rst) begin
        to_ifft_tlast_reg <= 0;
        mult_remain <= 0;
        state <= AWAIT_FRAME;
        mult_valid_proc_reg <= 0;
        frame_prep_done_corr <= 0;
        freq_to_phase_valid <= 0;
        ifft_data_after_coeff_to_sum <= 0;
        to_ifft_tvalid_reg <= 0;
        ready_ifft_reg <= 0;
        after_coeff_mult_done <= 0;
        phase_to_mult_reg <= 0;
        ifft_data_after_coeff <= 0;
        //ifft_data_after_coeff[1] <= 0;
        frame_prep_done_mwf <= 0;
        mwf_coef_req <= 0;
        m_corr_done <= 0;
        matrix_corr_do_div <= 0;
    //mwf_word_updated <= 0;
        mwf_matrix_w_updated <= 0;
        matrix_corr_est_mode <= 0;
        qL <= 0;
        mult_sel <= MULT_FREQ;
        qR <= 0;
        qL[CHANNEL_SIZE/2-1:0] <= {{INTEGER_BITS-1{1'b0}}, 1'b1, {DECIMAL_BITS{1'b0}}};
        qR[2*CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= {{INTEGER_BITS-1{1'b0}}, 1'b1, {DECIMAL_BITS{1'b0}}};
        words_aux <= 0;
        ready_reg <= 1;
        frame_vad <= 0;
        to_corr_ack <= 0;
        sample_addr <= 0;
        cont_frame <= 0;
        word_valid <= 0;
        phase_addr <= 0;
        frame_word_out <= 0;
        words_read <= 0;
    end else begin
        to_ifft_tlast_reg <= 0;
        to_ifft_tvalid_reg <= 0;
        mult_valid_proc_reg <= 0;
        frame_prep_done_corr <= 0;
        frame_prep_done_mwf <= 0;
        mwf_coef_req <= 0;
        m_corr_done <= 0;
        matrix_corr_do_div <= 0;
    //mwf_word_updated <= 0;
        mwf_matrix_w_updated <= 0;
        //freq_to_phase_valid <= 0;
        case(state)
            AWAIT_FRAME:  begin
                ready_reg <= 1;
                if (freq_wrd_valid_i) begin
                    //frame_word_in <= freq_word_i;
                    freq_word_to_phase[words_read*INPUT_WIDTH +: INPUT_WIDTH] <= freq_word_i;
                    freq_to_phase_valid[words_read*NUM_CHANNELS +: NUM_CHANNELS] <= {NUM_CHANNELS{1'b1}};
                    words_read <= words_read + 1;
                    frame_vad <= frame_vad | vad_i;
                    
                    phase_to_mult_reg <= phase_pattern_i;
                    if (words_read == CONC_MULTS-1 || sample_addr == NFFT-1) begin
                        ready_reg <= 0;
                        state <= PHASE_ADJUST;
                        finished_mults <= 0;
                    end else begin
                        sample_addr <= sample_addr + 1;
                    end
                    if (sample_addr < NFFT) begin
                        
                        if (phase_pattern_i == 0) begin
                            for (i=0; i<NUM_CHANNELS; i=i+1) begin
                                phase_to_mult_arr[words_read*INPUT_WIDTH+i*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_0_mem[sample_addr];
                            end
                            //phase_to_mult <= phase_0_mem[sample_addr];
                        end else if (phase_pattern_i == 1) begin
                            for (i=0; i<NUM_CHANNELS; i=i+1) begin
                                phase_to_mult_arr[words_read*INPUT_WIDTH+i*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_1_mem[sample_addr];
                            end
                            //phase_to_mult_arr[words_read*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_1_mem[sample_addr];
                            //phase_to_mult <= phase_1_mem[sample_addr];
                        end else if (phase_pattern_i == 2) begin
                            for (i=0; i<NUM_CHANNELS; i=i+1) begin
                                phase_to_mult_arr[words_read*INPUT_WIDTH+i*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_2_mem[sample_addr];
                            end
                            //phase_to_mult_arr[words_read*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_2_mem[sample_addr];
                            //phase_to_mult <= phase_2_mem[sample_addr];
                        end else begin
                            for (i=0; i<NUM_CHANNELS; i=i+1) begin
                                phase_to_mult_arr[words_read*INPUT_WIDTH+i*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_3_mem[sample_addr];
                            end
                            //phase_to_mult_arr[words_read*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_3_mem[sample_addr];
                            //phase_to_mult <= phase_3_mem[sample_addr];
                        end
                    end
                    //word_valid <= 1;
                end else begin
                    //word_valid <= 0;
                end
            end
            PHASE_ADJUST:  begin
                    for (i=0; i<CONC_MULTS*NUM_CHANNELS; i=i+1) begin
                        if (mult_ready[i] && freq_to_phase_valid[i]) begin
                            mult_in_a_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= freq_word_to_phase[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                            mult_in_b_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= phase_to_mult_arr[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                            mult_valid_proc_reg[i] <= 1;
                            freq_to_phase_valid[i] <= 0;
                        end
                        if (mult_done[i]) begin
                            finished_mults[i] <= 1;
                            stored_phase_res[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= mult_result[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                            //phase_addr
                        end
                    end
                    if (sample_addr < NFFT) begin
                        if (&(finished_mults)) begin
                            //finished_mults <= 0;
                            //frame_word_in <= ;
                            phase_addr_aux <= 0;
                            state <= STORE_PHASED;
                        end
                    end 
                
            end
            STORE_PHASED: begin
                freq_frame[phase_addr+phase_addr_aux] <= stored_phase_res[phase_addr_aux*CHANNEL_SIZE*NUM_CHANNELS +: CHANNEL_SIZE*NUM_CHANNELS];
                phase_addr_aux <= phase_addr_aux + 1;
                if (phase_addr_aux == words_read-1) begin
                    if (sample_addr == NFFT-1) begin
                        //phase_addr <= 0;
                        sample_addr <= 0;
                        cont_frame <= cont_frame + 1;
                        
                        if (matrix_corr_ready) begin
                            state <= AWAIT_CORR;
                            if (cont_frame == PRIOR_EST-1) begin  //trocar de volta para PRIOR_EST-1
                                matrix_corr_est_mode <= 1;
                                if (cont_frame == PRIOR_EST-1) begin
                                    matrix_corr_do_div <= 1;
                                end
                            end else begin
                                matrix_corr_est_mode <= 0;
                            end
                            mult_sel <= MULT_CORR;
                            frame_prep_done_corr <= 1;
                        end
                    end else begin
                        state <= AWAIT_FRAME;
                        sample_addr <= sample_addr + 1;
                        words_read <= 0;
                        ready_reg <= 1;
                        //phase_addr <= phase_addr + phase_addr_aux;
                        phase_addr <= sample_addr+1;
                    end
                end
            end
            AWAIT_CORR:  begin
                to_corr_ack <= 0;
                if (from_corr_req) begin
                    frame_word_out <= freq_frame[to_corr_address];
                    to_corr_ack <= 1;
                end
                if (matrix_corr_done) begin
                    m_corr_done <= 1;
                    if (matrix_corr_est_mode) begin
                        frame_prep_done_mwf <= 1;
                        mult_sel <= MULT_MWF;
                        state <= AWAIT_MWF;
                    end else begin
                        state <= AWAIT_qLR_MULT;
                        frame_word_out <= freq_frame[0];
                        mult_remain <= {(2*NUM_CHANNELS){1'b1}};
                        sample_addr <= 0;
                        ifft_data_after_coeff <= 0;
                        qLR_mult_done <= 0;
                        mult_sel <= MULT_FREQ;
                    end
                end
            end
            AWAIT_qLR_MULT: begin
                for (i=0; i<NUM_CHANNELS; i=i+1) begin
                        if (mult_ready[i] && mult_ready[i+NUM_CHANNELS] && mult_remain[i] && mult_remain[i+NUM_CHANNELS]) begin
                            mult_in_a_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, qL[i*CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]};
                            mult_in_a_proc_reg[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, qR[i*CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]};
                        
                            mult_in_b_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_word_out[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                            mult_in_b_proc_reg[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_word_out[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                        
                            mult_valid_proc_reg[i] <= 1;
                            mult_valid_proc_reg[i+NUM_CHANNELS] <= 1;
                            mult_remain[i] <= 0;
                            mult_remain[i+NUM_CHANNELS] <= 0;
                        end
                        if (mult_done[i]) begin
                            qLR_mult_done[i] <= 1;
                            //[i*CHANNEL_SIZE +: CHANNEL_SIZE/2]
                            ifft_data_after_coeff_to_sum[i*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= mult_result[i*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            ifft_data_after_coeff_to_sum[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= mult_result[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                            //ifft_data_after_coeff[CHANNEL_SIZE/2-1:0] <= ifft_data_after_coeff[CHANNEL_SIZE/2-1:0] + mult_result[i*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            //ifft_data_after_coeff[CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mult_result[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end
                        if (mult_done[i+NUM_CHANNELS]) begin
                            qLR_mult_done[i+NUM_CHANNELS] <= 1;
                            ifft_data_after_coeff_to_sum[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            ifft_data_after_coeff_to_sum[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                            //ifft_data_after_coeff[CHANNEL_SIZE +: CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE +: CHANNEL_SIZE/2] + mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            //ifft_data_after_coeff[CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end
                end
                if (&qLR_mult_done) begin
                        mult_remain[1:0] <= 2'b11;
                        ifft_data_after_coeff <= ifft_data_after_coeff_sum;
                        //ifft_data_after_coeff_sum
                        state <= TO_IFFT_CONJ_PHASE_MULT;
                        if (phase_to_mult_reg == 0) begin
                            phase_to_mult <= phase_0_mem[sample_addr];
                        end else if (phase_to_mult_reg == 1) begin
                            phase_to_mult <= phase_1_mem[sample_addr];
                        end else if (phase_to_mult_reg == 2) begin
                            phase_to_mult <= phase_2_mem[sample_addr];
                        end else begin
                            phase_to_mult <= phase_3_mem[sample_addr];
                        end
                end
            end
            AWAIT_MWF:  begin
                if (mwf_done) begin
                    mwf_coef_req <= 1;
                    state <= AWAIT_CONJ_COEF_MULT;
                    mult_remain[2*NUM_CHANNELS-1:0] <= {(2*NUM_CHANNELS){1'b1}};
                    ifft_data_after_coeff <= 0;
                    //ifft_data_after_coeff[1] <= 0;
                    mult_sel <= MULT_FREQ;
                    after_coeff_mult_done <= 0;
                    frame_word_out <= freq_frame[0];
                    sample_addr <= 0;
                end
            end
            AWAIT_CONJ_COEF_MULT:  begin
                if (mwf_coef_ack) begin
                    for (i=0; i<NUM_CHANNELS; i=i+1) begin
                        if (mult_ready[i] && mult_ready[i+NUM_CHANNELS] && mult_remain[i] && mult_remain[i+NUM_CHANNELS]) begin
                            mult_in_a_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= {~coef_sample[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]+1, coef_sample[i*CHANNEL_SIZE +: CHANNEL_SIZE/2]};
                            mult_in_a_proc_reg[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE] <= {~coef_sample[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]+1, coef_sample[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2]};
                        
                            mult_in_b_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_word_out[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                            mult_in_b_proc_reg[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_word_out[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                        
                            mult_valid_proc_reg[i] <= 1;
                            mult_valid_proc_reg[i+NUM_CHANNELS] <= 1;
                            mult_remain[i+NUM_CHANNELS] <= 0;
                            mult_remain[i] <= 0;
                        end
                    end
                end
                    for (i=0; i<NUM_CHANNELS; i=i+1) begin
                        if (mult_done[i]) begin
                            //fazer soma combinacional aqui!
                            after_coeff_mult_done[i] <= 1;
                            ifft_data_after_coeff[CHANNEL_SIZE/2-1:0] <= ifft_data_after_coeff[CHANNEL_SIZE/2-1:0] + mult_result[i*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            ifft_data_after_coeff[CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mult_result[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end
                        if (mult_done[i+NUM_CHANNELS]) begin
                            
                            after_coeff_mult_done[i+NUM_CHANNELS] <= 1;
                            ifft_data_after_coeff[CHANNEL_SIZE +: CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE +: CHANNEL_SIZE/2] + mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            ifft_data_after_coeff[CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= ifft_data_after_coeff[CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + mult_result[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end
                    end
                    if (&after_coeff_mult_done) begin
                        state <= TO_IFFT_CONJ_PHASE_MULT;
                        mult_remain[1:0] <= 2'b11;
                        if (phase_to_mult_reg == 0) begin
                            phase_to_mult <= phase_0_mem[sample_addr];
                        end else if (phase_to_mult_reg == 1) begin
                            phase_to_mult <= phase_1_mem[sample_addr];
                        end else if (phase_to_mult_reg == 2) begin
                            phase_to_mult <= phase_2_mem[sample_addr];
                        end else begin
                            phase_to_mult <= phase_3_mem[sample_addr];
                        end
                    end
            end
            TO_IFFT_CONJ_PHASE_MULT: begin
                for (i=0; i<2; i=i+1) begin
                    if (mult_ready[i] && mult_remain[i] == 1) begin
                        mult_in_a_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= ifft_data_after_coeff[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                        mult_in_b_proc_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= {~phase_to_mult[CHANNEL_SIZE-1:CHANNEL_SIZE/2]+1, phase_to_mult[CHANNEL_SIZE/2-1:0]};
                        mult_valid_proc_reg[i] <= 1;
                        mult_remain[i] <= 0;
                    end
                    if (mult_done[i]) begin
                        after_coeff_mult_done[i] <= 1;
                        to_ifft_tdata_reg[i*CHANNEL_SIZE +: CHANNEL_SIZE] <= mult_result[i*CHANNEL_SIZE +: CHANNEL_SIZE];
                    end
                end
                if (&(after_coeff_mult_done[1:0])) begin
                    to_ifft_mem[sample_addr] <= to_ifft_tdata_reg;
                    if (sample_addr == NFFT-1) begin
                        sample_addr <= 0;
                        state <= READY_FOR_IFFT;
                        ready_ifft_reg <= 1;
                    end else begin
                        sample_addr <= sample_addr + 1;
                        after_coeff_mult_done <= 0;
                        ifft_data_after_coeff <= 0;
                        frame_word_out <= freq_frame[sample_addr+1];
                        if (cont_frame <= PRIOR_EST-1) begin
                            qLR_mult_done <= 0;
                            state <= AWAIT_qLR_MULT;
                            mult_remain <= {(2*NUM_CHANNELS){1'b1}};
                        end else begin
                            mwf_coef_req <= 1;
                            state <= AWAIT_CONJ_COEF_MULT;
                            mult_remain[2*NUM_CHANNELS-1:0] <= {(2*NUM_CHANNELS){1'b1}};
                        end
                        
                    end
                end
            end
            READY_FOR_IFFT: begin
                if (ifft_core_ready_i) begin
                    to_ifft_tdata_reg <= to_ifft_mem[sample_addr];
                    to_ifft_tvalid_reg <= 1;
                    if (sample_addr == NFFT-1) begin
                        to_ifft_tlast_reg <= 1;
                        sample_addr <= 0;
                        state <= AWAIT_FRAME;
                        ready_ifft_reg <= 0;
                        words_read <= 0;
                        ready_reg <= 1;
                    end else begin
                        sample_addr <= sample_addr + 1;
                        to_ifft_tlast_reg <= 0;
                    end
                end
            end
            default:
                state <= AWAIT_FRAME;
        endcase
    end
   end
   
   always @(state, qLR_mult_done/*, after_coeff_mult_done*/) begin
    if (state == AWAIT_qLR_MULT) begin
        if (&qLR_mult_done) begin
            ifft_data_after_coeff_sum = 0;
            for (i=0; i<NUM_CHANNELS; i=i+1) begin
                ifft_data_after_coeff_sum[CHANNEL_SIZE/2-1:0] = ifft_data_after_coeff_sum[CHANNEL_SIZE/2-1:0] + ifft_data_after_coeff_to_sum[i*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                ifft_data_after_coeff_sum[CHANNEL_SIZE-1:CHANNEL_SIZE/2] = ifft_data_after_coeff_sum[CHANNEL_SIZE-1:CHANNEL_SIZE/2] + ifft_data_after_coeff_to_sum[i*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                
                ifft_data_after_coeff_sum[CHANNEL_SIZE+CHANNEL_SIZE/2-1:CHANNEL_SIZE] = ifft_data_after_coeff_sum[CHANNEL_SIZE+CHANNEL_SIZE/2-1:CHANNEL_SIZE] + ifft_data_after_coeff_to_sum[(i+NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                ifft_data_after_coeff_sum[CHANNEL_SIZE+CHANNEL_SIZE-1:CHANNEL_SIZE+CHANNEL_SIZE/2] = ifft_data_after_coeff_sum[CHANNEL_SIZE-1:CHANNEL_SIZE/2] + ifft_data_after_coeff_to_sum[(i+NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
            end
        end
    end
//    if (state == TO_IFFT_CONJ_PHASE_MULT) begin
//        if (&(after_coeff_mult_done[1:0]) == 2'b11) begin
            
//        end
//    end
   end
   
   always @(mult_sel, mult_valid_corr, mult_valid_proc, mult_valid_mwf) begin
        case(mult_sel)
            MULT_FREQ: begin
                mult_valid_selected = mult_valid_proc;
            end
            MULT_CORR: begin
                mult_valid_selected = mult_valid_corr;
            end
            MULT_MWF: begin
                mult_valid_selected = mult_valid_mwf;
            end 
            default:
                mult_valid_selected = 0;
        endcase
   end
   always @(mult_sel, mult_in_a_corr, mult_in_a_proc, mult_in_a_mwf) begin
        case(mult_sel)
            MULT_FREQ: begin
                mult_in_a_selected = mult_in_a_proc;
            end
            MULT_CORR: begin
                mult_in_a_selected = mult_in_a_corr;
            end
            MULT_MWF: begin
                mult_in_a_selected = mult_in_a_mwf;
            end 
            default:
                mult_in_a_selected = 0;
        endcase
   end
   always @(mult_sel, mult_in_b_corr, mult_in_b_proc, mult_in_b_mwf) begin
        case(mult_sel)
            MULT_FREQ: begin
                mult_in_b_selected = mult_in_b_proc;
            end
            MULT_CORR: begin
                mult_in_b_selected = mult_in_b_corr;
            end
            MULT_MWF: begin
                mult_in_b_selected = mult_in_b_mwf;
            end 
            default:
                mult_in_b_selected = 0;
        endcase
   end
               
   matrix_corr_alt #(
    .NFFT(NFFT),
    .INTEGER_BITS(INTEGER_BITS),
    .NUM_PATTERNS(NUM_PATTERNS),
    .DECIMAL_BITS(DECIMAL_BITS),
    .CHANNEL_SIZE(CHANNEL_SIZE),
    .NUM_CHANNELS(NUM_CHANNELS),
    .MULT_AMOUNT(MULT_AMOUNT),
    .PRIOR_EST(PRIOR_EST),
    .INPUT_WIDTH(INPUT_WIDTH),
    .MATRIX_WORD_SIZE(CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS)
  ) matrix_inst (
    .clk(clk),
    .rst(rst),
    .sample_i(frame_word_out),
    .valid_frame_i(frame_prep_done_corr),
    .sample_ack_i(to_corr_ack),
    .words_requested_i(mwf_w_update_req),
    .word_req_i(to_mwf_address),
    .word_ack_o(mwf_word_updated),
    .req_rx_word_o(corr_to_mwf_rx_word),
    .req_ry_word_o(corr_to_mwf_ry_word),
    .req_new_sample_o(from_corr_req),
    .sample_req_o(to_corr_address),
    .vad_i(frame_vad),
    .estimation_mode_i(matrix_corr_est_mode),
    .do_division_i(matrix_corr_do_div),
    .check_empty(check_empty),
    .corr_mult_done(mult_done),
    .corr_mult_ready(mult_ready),
    .corr_mult_result(mult_result),
    .corr_to_mult_valid(mult_valid_corr),
    .corr_to_mult_real_mult(),
    .corr_to_mult_in_a(mult_in_a_corr),
    .corr_to_mult_in_b(mult_in_b_corr),
    .div_done_o(matrix_corr_div_done),
    .ready_o(matrix_corr_ready),
    .done_o(matrix_corr_done)
    );
    
    mwf_alt # (
    .NFFT(NFFT),
    .NUM_CHANNELS(NUM_CHANNELS), 
    .CHANNEL_SIZE(CHANNEL_SIZE),
    .INPUT_WIDTH(INPUT_WIDTH),
    .MULT_AMOUNT(MULT_AMOUNT),
    .DECIMAL_BITS(DECIMAL_BITS)
) mwf_inst (
    .clk(clk),
    .rst(rst),
    .frame_valid_i(frame_prep_done_mwf),
    .rv_packed_word(0),
    .ry_packed_word(corr_to_mwf_ry_word),
    .rx_packed_word(corr_to_mwf_rx_word),
    .word_updated_i(mwf_word_updated),
    .sampe_addr_o(to_mwf_address),
    .update_word_o(mwf_w_update_req),
    .matrices_updated_i(m_corr_done),
    .qL_packed(qL),
    .qR_packed(qR),
    .mwf_mult_done(mult_done),
    .mwf_mult_ready(mult_ready),
    .mwf_mult_result(mult_result),
    .mwf_to_mult_valid(mult_valid_mwf),
    .mwf_to_mult_real_mult(),
    .mwf_to_mult_in_a(mult_in_a_mwf),
    .mwf_to_mult_in_b(mult_in_b_mwf),
    .coeff_packed(coef_sample),
    .coeff_sample_request_i(mwf_coef_req),
    .coeff_sample_addr_i(sample_addr),
    .coeff_sample_ack_o(mwf_coef_ack),
    .ready_o(mwf_ready),
    .done_o(mwf_done)
 );
endmodule
