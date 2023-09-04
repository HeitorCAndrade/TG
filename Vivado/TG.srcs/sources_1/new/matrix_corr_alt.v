`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 22.08.2023 15:38:11
// Design Name: 
// Module Name: matrix_corr_alt
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


module matrix_corr_alt #(
    parameter NFFT = 256,
    parameter INTEGER_BITS = 12,
    parameter NUM_PATTERNS = 4,
    parameter DECIMAL_BITS = 20,
    parameter CHANNEL_SIZE = 64,
    parameter NUM_CHANNELS = 6,
    parameter MULT_AMOUNT = 32,
    parameter PRIOR_EST = 1001,
    parameter INPUT_WIDTH = CHANNEL_SIZE*NUM_CHANNELS,
    parameter MATRIX_WORD_SIZE = CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS
)(
    input clk,
    input rst,
    input wire signed [INPUT_WIDTH-1:0] sample_i,
    input wire valid_frame_i,
    input wire sample_ack_i,
    input wire words_requested_i,
    input wire [7:0] word_req_i,
    output wire word_ack_o,
    output wire [MATRIX_WORD_SIZE-1:0] req_rx_word_o,
    output wire [MATRIX_WORD_SIZE-1:0] req_ry_word_o,
    output wire req_new_sample_o,
    output wire [7:0] sample_req_o,
    input wire vad_i,
    input wire estimation_mode_i,
    input wire do_division_i,
    input wire check_empty,
    input wire [MULT_AMOUNT-1:0] corr_mult_done,
    input wire [MULT_AMOUNT-1:0] corr_mult_ready,
    input wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] corr_mult_result,
    output wire [MULT_AMOUNT-1:0] corr_to_mult_valid,
    output wire [MULT_AMOUNT-1:0] corr_to_mult_real_mult,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] corr_to_mult_in_a,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] corr_to_mult_in_b,
    output wire div_done_o,
//    output wire signed [2*NFFT*INPUT_WIDTH-1:0]compressed_rv_o,
//    output wire rv_valid_o,
//    output wire signed [2*NFFT*INPUT_WIDTH-1:0]compressed_ry_o,
//    output wire ry_valid_o,
//    output wire signed [2*NFFT*INPUT_WIDTH-1:0]compressed_rx_o,
//    output wire rx_valid_o,
    output wire ready_o,
    output wire done_o
    );
    
    genvar i;
    integer channel, channel2;
    localparam READY = 0,
               FIRST_MULT = 1,
               EST_0_SUM = 2,
               EST_1_MULT_1 = 3,
               EST_1_MULT_2 = 4,
               UPDATE_X = 5,
               DO_DIV = 6,
               FETCH_SAMPLE = 7,
               SEND_REQ_WORDS = 8,
               SKIP_MATRIX = 9;
               
//localparam CONCURRENT_SAMPLES = 5;
localparam CONCURRENT_MULTS = MULT_AMOUNT-2;
//localparam MATRIX_WORD_SIZE = CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS;
localparam DIV_AMOUNT = 12;
assign corr_to_mult_real_mult = 0;       

reg signed [CHANNEL_SIZE/2-1:0] smooth_coef = 32'h000FCA18;
reg signed [CHANNEL_SIZE/2-1:0] smooth_coef_compl = 32'h000035E7; //(1-smooth_coef)

reg [3:0] state;
reg ready_reg; 
reg done_reg;
reg noise_only;
reg est_mode_reg;
reg do_div_reg;
reg check_empty_reg;


reg req_new_sample_reg;
reg [7:0] sample_counter;

reg [CHANNEL_SIZE-1:0] frame_sample_reg [0:NUM_CHANNELS-1];
wire [CHANNEL_SIZE-1:0] frame_y_inv [0:NUM_CHANNELS-1];
wire [CHANNEL_SIZE-1:0] frame_y_conj [0:NUM_CHANNELS-1];

reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] corr_to_mult_in_a_reg;
reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] corr_to_mult_in_b_reg;
reg [MULT_AMOUNT-1:0] corr_to_mult_valid_reg;
reg [2:0] channel_1_reg, channel_2_reg;
reg [6:0] mult_completed;
reg [MULT_AMOUNT-1:0] finished_mults;
reg [NUM_CHANNELS-1:0] reamining_mult;
reg [CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS-1:0] frame_y_new;

reg channel_addr_reg;
reg [MATRIX_WORD_SIZE-1:0] matrix_word, matrix_word_2;
reg [4:0] matrix_word_offset;
reg word_updated;
reg [CHANNEL_SIZE/2-1:0] cont_ry;
reg [CHANNEL_SIZE/2-1:0] cont_rv;
reg [CHANNEL_SIZE/4-1:0] cont_rx;
(* ram_style = "block" *) reg [MATRIX_WORD_SIZE-1:0] rv [0:NFFT-1];
(* ram_style = "block" *) reg [MATRIX_WORD_SIZE-1:0] ry [0:NFFT-1];
(* ram_style = "block" *) reg [MATRIX_WORD_SIZE-1:0] rx [0:NFFT-1];

reg [CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS-1:0] after_est_mult;

reg [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_in_a_re;
reg [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_in_b_re;
reg [DIV_AMOUNT-1:0] div_in_start_re;
wire [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_out_val_re;
wire [DIV_AMOUNT-1:0] div_out_busy_re;
wire [DIV_AMOUNT-1:0] div_out_done_re;
wire [DIV_AMOUNT-1:0] div_out_valid_re;
wire [DIV_AMOUNT-1:0] div_out_dbz_re;
wire [DIV_AMOUNT-1:0] div_out_ovf_re;

reg [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_in_a_im;
reg [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_in_b_im;
reg [DIV_AMOUNT-1:0] div_in_start_im;
wire [DIV_AMOUNT*CHANNEL_SIZE/2-1:0] div_out_val_im;
wire [DIV_AMOUNT-1:0] div_out_busy_im;
wire [DIV_AMOUNT-1:0] div_out_done_im;
wire [DIV_AMOUNT-1:0] div_out_valid_im;
wire [DIV_AMOUNT-1:0] div_out_dbz_im;
wire [DIV_AMOUNT-1:0] div_out_ovf_im;

reg word_ack_reg;
reg div_done_reg;



generate
    for (i=0; i<NUM_CHANNELS; i=i+1) begin
        assign frame_y_inv[i] = frame_sample_reg[i];
        assign frame_y_conj[i][CHANNEL_SIZE/2-1:0] = frame_sample_reg[i][CHANNEL_SIZE/2-1:0];
        assign frame_y_conj[i][CHANNEL_SIZE-1:CHANNEL_SIZE/2] = frame_sample_reg[i][CHANNEL_SIZE-1:CHANNEL_SIZE/2];
    end
    for (i=0; i<DIV_AMOUNT; i=i+1) begin
        div # (
            .WIDTH(CHANNEL_SIZE/2),
            .FBITS(DECIMAL_BITS)
        ) div_re (
            .clk(clk),
            .rst(rst),
            .start(div_in_start_re[i]),
            .busy(div_out_busy_re[i]),
            .done(div_out_done_re[i]),
            .valid(div_out_valid_re[i]),
            .dbz(div_out_dbz_re[i]),
            .ovf(div_out_ovf_re[i]),
            .a(div_in_a_re[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2]),
            .b(div_in_b_re[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2]),
            .val(div_out_val_re[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2])
         );       
         div # (
            .WIDTH(CHANNEL_SIZE/2),
            .FBITS(DECIMAL_BITS)
        ) div_im (
            .clk(clk),
            .rst(rst),
            .start(div_in_start_im[i]),
            .busy(div_out_busy_im[i]),
            .done(div_out_done_im[i]),
            .valid(div_out_valid_im[i]),
            .dbz(div_out_dbz_im[i]),
            .ovf(div_out_ovf_im[i]),
            .a(div_in_a_im[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2]),
            .b(div_in_b_im[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2]),
            .val(div_out_val_im[i*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2])
         );
    end
endgenerate
assign done_o = done_reg;
assign ready_o = ready_reg;
assign sample_req_o = sample_counter;
assign req_rx_word_o = matrix_word;
assign req_ry_word_o = matrix_word_2;
assign word_ack_o = word_ack_reg;
assign req_new_sample_o = req_new_sample_reg;
assign div_done_o = div_done_reg;

assign corr_to_mult_valid = corr_to_mult_valid_reg;
assign corr_to_mult_in_a = corr_to_mult_in_a_reg;
assign corr_to_mult_in_b = corr_to_mult_in_b_reg;
always @(posedge(clk)) begin
    ready_reg <= 0;
    done_reg <= 0;
    div_in_start_re <= 0;
    div_in_start_im <= 0;
    corr_to_mult_valid_reg <= 0;
    req_new_sample_reg <= 0;
    word_ack_reg <= 0;
    if (rst) begin
        state <= READY;
        check_empty_reg <= 0;
        div_done_reg <= 0;
        do_div_reg <= 0;
        sample_counter <= 0;
        cont_rv <= 0;
        cont_ry <= 0;
        cont_rx <= 0;
        matrix_word <= 0;
        matrix_word_2 <= 0;
    end else begin 
        case (state)
            READY: begin
                ready_reg <= 1;
                req_new_sample_reg <= 0;
                if (valid_frame_i) begin
                    if (check_empty) begin
                        est_mode_reg <= estimation_mode_i;
                        do_div_reg <= do_division_i;
                        noise_only <= ~vad_i;
                        req_new_sample_reg <= 1;
                        //check_empty_reg <= check_empty;
                        sample_counter <= 0;
                        state <= FETCH_SAMPLE;
                    end else begin
                        done_reg <= 1;
                    end
                end else if (words_requested_i) begin
                    sample_counter <= word_req_i;
                    state <= SEND_REQ_WORDS;
                end
            end
            SKIP_MATRIX: begin
            end
            SEND_REQ_WORDS: begin
                matrix_word <= rx[sample_counter];
                matrix_word_2 <= ry[sample_counter];
                word_ack_reg <= 1;
                state <= READY;
            end
            FETCH_SAMPLE: begin
                if (sample_ack_i) begin
                    for (channel=0; channel<NUM_CHANNELS; channel=channel+1) begin
                        frame_sample_reg[channel] <= sample_i[channel*CHANNEL_SIZE +: CHANNEL_SIZE];
                    end
                    channel_1_reg <= 0;
                    mult_completed <= 0;
                    finished_mults <= 0;
                    matrix_word_offset <= 0;
                    reamining_mult <= {NUM_CHANNELS{1'b1}};
                    state <= FIRST_MULT;
                end
            end
            FIRST_MULT: begin
                for (channel=0; channel<NUM_CHANNELS; channel=channel+1) begin
                    if (corr_mult_ready[channel] && reamining_mult[channel]) begin
                        corr_to_mult_in_a_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_y_inv[channel_1_reg];
                        corr_to_mult_in_b_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_y_conj[channel];
                        corr_to_mult_valid_reg[channel] <= 1;
                        reamining_mult[channel] <= 0;
                    end
                    if (corr_mult_done[channel]) begin
                        finished_mults[channel] <= 1;
                        frame_y_new[(matrix_word_offset+channel)*CHANNEL_SIZE +: CHANNEL_SIZE] <= corr_mult_result[channel*CHANNEL_SIZE +: CHANNEL_SIZE];
                    end
                end
                if (channel_1_reg < 5) begin
                    if (&(finished_mults[NUM_CHANNELS-1:0])) begin
                        channel_1_reg <= channel_1_reg + 1;
                        matrix_word_offset <= matrix_word_offset + NUM_CHANNELS;
                        finished_mults <= 0;
                        reamining_mult <= {NUM_CHANNELS{1'b1}};
                    end
                end else begin
                    if (&(finished_mults[NUM_CHANNELS-1:0])) begin
                        if (noise_only) begin
                            matrix_word <= rv[sample_counter];
                        end else begin
                            matrix_word <= ry[sample_counter];
                        end
                        if (est_mode_reg == 0) begin 
                            channel_addr_reg <= 0;
                            word_updated <= 0;
                            state <= EST_0_SUM;
                        end else begin
                            state <= EST_1_MULT_1;
                            matrix_word_offset <= 0;
                            finished_mults <= 0;
                            channel_1_reg <= 0;
                        end
                    end
                end
            end
            EST_0_SUM: begin
                word_updated <= 1;
                for (channel=0; channel<NUM_CHANNELS; channel=channel+1) begin
                    for (channel2=0; channel2<NUM_CHANNELS; channel2=channel2+1) begin
                        matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE +: CHANNEL_SIZE/2] + frame_y_new[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                        matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + frame_y_new[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                    end
                end
                if (word_updated) begin
                    if (noise_only) begin
                        rv[sample_counter] <= matrix_word;
                    end else begin
                        ry[sample_counter] <= matrix_word;
                    end
                    if (sample_counter == NFFT-1) begin
                        if (do_div_reg) begin
                            sample_counter <= 0;
                            matrix_word <= rv[0];
                            matrix_word_2 <= ry[0];
                            finished_mults <= 0;
                            word_updated <= 1;
                            state <= DO_DIV;
                        end else begin
                            state <= READY;
                            done_reg <= 1;
                        end
                        if (noise_only) begin
                            cont_rv <= cont_rv + 1;
                        end else begin
                            cont_ry <= cont_ry + 1;
                        end
                    end else begin 
                        sample_counter <= sample_counter + 1;
                        req_new_sample_reg <= 1;
                        state <= FETCH_SAMPLE;
                    end
                end
            end
            EST_1_MULT_1: begin
                for (channel=0; channel<CONCURRENT_MULTS; channel=channel+1) begin
                    if (corr_mult_ready[channel]) begin
                        corr_to_mult_valid_reg[channel] <= 1; //smooth_coef_compl
                        corr_to_mult_in_b_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, smooth_coef_compl};
                        if (channel_1_reg < 5 || channel < NUM_CHANNELS) begin
                            corr_to_mult_in_a_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= frame_y_new[(matrix_word_offset+channel)*CHANNEL_SIZE +: CHANNEL_SIZE];
                        end else begin
                            corr_to_mult_in_a_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= 0;
                        end
                    end
                    if (corr_mult_done[channel] && matrix_word_offset+channel < NUM_CHANNELS*NUM_CHANNELS) begin
                        finished_mults[channel] <= 1;
                        after_est_mult[(matrix_word_offset+channel)*CHANNEL_SIZE +: CHANNEL_SIZE] <= corr_mult_result[channel*CHANNEL_SIZE +: CHANNEL_SIZE];
                    end
                end 
                if (channel_1_reg < 5) begin
                    if (&(finished_mults[CONCURRENT_MULTS-1:0])) begin
                        channel_1_reg <= channel_1_reg + 5;
                        matrix_word_offset <= 30;
                        finished_mults <= 0;
                    end
                end else begin
                    if (&(finished_mults[NUM_CHANNELS-1:0])) begin
                        state <= EST_1_MULT_2;
                        channel_1_reg <= 0;
                        matrix_word_offset <= 0;
                        finished_mults <= 0;
                    end
                end
            end
            EST_1_MULT_2: begin
                for (channel=0; channel<CONCURRENT_MULTS; channel=channel+1) begin
                    if (corr_mult_ready[channel]) begin
                        corr_to_mult_valid_reg[channel] <= 1;
                        corr_to_mult_in_b_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, smooth_coef};
                        if (channel_1_reg < 5 || channel < NUM_CHANNELS) begin
                            corr_to_mult_in_a_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= matrix_word[(matrix_word_offset+channel)*CHANNEL_SIZE +: CHANNEL_SIZE];
                        end else begin
                            corr_to_mult_in_a_reg[channel*CHANNEL_SIZE +: CHANNEL_SIZE] <= 0;
                        end
                    end
                    if (corr_mult_done[channel] && matrix_word_offset+channel < NUM_CHANNELS*NUM_CHANNELS) begin//after_est_mult
                        finished_mults[channel] <= 1;
                        matrix_word[(matrix_word_offset+channel)*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= corr_mult_result[channel*CHANNEL_SIZE +: CHANNEL_SIZE/2] + after_est_mult[channel*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                        matrix_word[(matrix_word_offset+channel)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= corr_mult_result[channel*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + after_est_mult[channel*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                    end
                end 
                if (channel_1_reg < 5) begin
                    if (&(finished_mults[CONCURRENT_MULTS-1:0])) begin
                        channel_1_reg <= channel_1_reg + 5;
                        matrix_word_offset <= 30;
                        finished_mults <= 0;
                    end
                end else begin
                    if (&(finished_mults[NUM_CHANNELS-1:0])) begin
                        if (noise_only) begin
                            rv[sample_counter] <= matrix_word;
                        end else begin
                            ry[sample_counter] <= matrix_word;
                        end
                        matrix_word <= rv[sample_counter];
                        matrix_word_2 <= ry[sample_counter];
                        word_updated <= 0;
                        state <= UPDATE_X;
                        finished_mults <= 0;
                    end
                end
            end
            UPDATE_X: begin
                word_updated <= 1;
                for (channel=0; channel<NUM_CHANNELS; channel=channel+1) begin
                    for (channel2=0; channel2<NUM_CHANNELS; channel2=channel2+1) begin
                        matrix_word_2[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= $signed(matrix_word_2[channel*NUM_CHANNELS*CHANNEL_SIZE +: CHANNEL_SIZE/2]) - $signed(matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE +: CHANNEL_SIZE/2]);
                        matrix_word_2[channel*NUM_CHANNELS*CHANNEL_SIZE+channel2*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= $signed(matrix_word_2[channel*NUM_CHANNELS*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]) - $signed(matrix_word[channel*NUM_CHANNELS*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2]);
                    end
                end
                if (word_updated) begin
                    rx[sample_counter] <= matrix_word;
                    if (sample_counter == NFFT-1) begin
                        state <= READY;
                        done_reg <= 1;
                    end else begin 
                        sample_counter <= sample_counter + 1;
                        req_new_sample_reg <= 1;
                        state <= FETCH_SAMPLE;
                    end
                end
            end
            DO_DIV: begin
                if (word_updated) begin 
                    for (channel=0; channel < DIV_AMOUNT; channel=channel+1) begin
                        if (div_out_busy_re[channel]==0 && div_out_busy_im[channel]==0) begin
                            div_in_start_re[channel] <= 1;
                            if (channel < NUM_CHANNELS) begin
                                div_in_a_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= matrix_word[channel*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                                div_in_a_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= matrix_word[channel*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                                div_in_b_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= cont_rv;
                                div_in_b_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= cont_rv;
                            end else begin 
                                div_in_a_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= matrix_word_2[(channel-NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                                div_in_a_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= matrix_word_2[(channel-NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                                div_in_b_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= cont_ry;
                                div_in_b_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2] <= cont_ry;
                            end
                        end
                        if (div_out_done_re[channel] && div_out_valid_re[channel]) begin
                            finished_mults[channel] <= 1;
                            if (channel < NUM_CHANNELS) begin
                                matrix_word[channel*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= div_out_val_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2];
                            end else begin
                                matrix_word_2[(channel-NUM_CHANNELS)*CHANNEL_SIZE +: CHANNEL_SIZE/2] <= div_out_val_re[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2];
                            end
                        end
                        if (div_out_done_im[channel] && div_out_valid_im[channel]) begin
                            finished_mults[channel+DIV_AMOUNT] <= 1;
                            if (channel < NUM_CHANNELS) begin
                                matrix_word[channel*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= div_out_val_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2];
                            end else begin
                                matrix_word_2[(channel-NUM_CHANNELS)*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] <= div_out_val_im[channel*(CHANNEL_SIZE/2) +: CHANNEL_SIZE/2];
                            end
                        end
                    end
                end else begin
                    matrix_word <= rv[sample_counter];
                    matrix_word_2 <= ry[sample_counter];
                    word_updated <= 1;
                end
                if (&(finished_mults[2*DIV_AMOUNT-1:0])) begin
                    finished_mults <= 0;
                    word_updated <= 0;
                    rv[sample_counter] <= matrix_word;
                    ry[sample_counter] <= matrix_word_2;
                    sample_counter <= sample_counter + 1;
                    if (sample_counter == NFFT-1) begin
                        done_reg <= 1;
                        div_done_reg <= 1;
                        state <= READY;
                    end 
                end
            end
            default:
                state <= READY;
        endcase
    end
end
    
endmodule
