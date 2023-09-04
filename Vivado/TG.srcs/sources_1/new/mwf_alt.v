`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 21.08.2023 22:26:00
// Design Name: 
// Module Name: mwf_alt
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


module mwf_alt # (
    parameter NFFT = 256,   
    parameter Fs = 16000,   
    parameter NUM_CHANNELS = 6, 
    parameter CHANNEL_SIZE = 64,
    parameter INPUT_WIDTH = CHANNEL_SIZE*NUM_CHANNELS,
    parameter MULT_AMOUNT = 32,
    parameter DECIMAL_BITS = 20
)
(
    input clk,
    input rst,
    //input wire [NFFT*INPUT_WIDTH-1:0] frame_packed,
    input wire frame_valid_i,
//    input wire [NUM_CHANNELS*NFFT*INPUT_WIDTH-1:0] rv_packed,
//    input wire [NUM_CHANNELS*NFFT*INPUT_WIDTH-1:0] ry_packed,
//    input wire [NUM_CHANNELS*NFFT*INPUT_WIDTH-1:0] rx_packed,
    input wire [NUM_CHANNELS*INPUT_WIDTH-1:0] rv_packed_word,
    input wire [NUM_CHANNELS*INPUT_WIDTH-1:0] ry_packed_word,
    input wire [NUM_CHANNELS*INPUT_WIDTH-1:0] rx_packed_word,
    input wire word_updated_i,
    output wire sampe_addr_o,
    output wire update_word_o,
    input wire matrices_updated_i,
    input wire signed [INPUT_WIDTH/2-1:0] qL_packed,
    input wire signed [INPUT_WIDTH/2-1:0] qR_packed,
    input wire [MULT_AMOUNT-1:0] mwf_mult_done,
    input wire [MULT_AMOUNT-1:0] mwf_mult_ready,
    input wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mwf_mult_result,
    output wire [MULT_AMOUNT-1:0] mwf_to_mult_valid,
    output wire [MULT_AMOUNT-1:0] mwf_to_mult_real_mult,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mwf_to_mult_in_a,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] mwf_to_mult_in_b,
    output wire [INPUT_WIDTH*2-1:0] coeff_packed,
    input wire coeff_sample_request_i,
    input wire [7:0] coeff_sample_addr_i,
    output wire coeff_sample_ack_o,
    output wire ready_o,
    output wire done_o
 );
 genvar i, j;
 integer channel_1, channel_2;
 localparam SAMPLE_USE = (NFFT/2)+1;
 localparam REDUCED_MULT_AMOUNT = MULT_AMOUNT-8; //para facilitar controle das somas, mudar para "MULT_AMOUNT" em eventual otimização
 
 //states
 localparam READY = 0,
            R_PREP = 1,
            V1_MULT_1 = 2,
            V1_MULT_2 = 3,
            V1_SUM = 4,
            M1_MULT = 5,
            MWF_CALC = 6,
            COEFF_UPDATE = 7,
            FETCH_M_WORDS = 8,
            SEND_COEFF_SAMPLE = 9;
 
 assign mwf_to_mult_real_mult = 0;
 reg [CHANNEL_SIZE/2-1:0] BETA = 32'h000051EB;
 reg [CHANNEL_SIZE/2-1:0] GAMMA = 0;
 
 reg [7:0] sample_counter;
 reg [3:0] state;
 reg ready_reg;
 reg frame_valid_reg;
 reg done_reg;
 
 wire [CHANNEL_SIZE-1:0] ry_w [0:NUM_CHANNELS*NUM_CHANNELS-1];
 wire [CHANNEL_SIZE-1:0] rx_w [0:NUM_CHANNELS*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS-1:0] matrix_word_rv;
 reg [CHANNEL_SIZE*NUM_CHANNELS*NUM_CHANNELS-1:0] matrix_word_ry;
 
 reg [CHANNEL_SIZE-1:0] rxx [0:2*NUM_CHANNELS-1][0:2*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE-1:0] ryy [0:2*NUM_CHANNELS-1][0:2*NUM_CHANNELS-1];
 
 wire [INPUT_WIDTH/2-1:0] qLR [0:2*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE-1:0] qLR_BETA_mult_result [0:2*NUM_CHANNELS-1];
 reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] mwf_to_mult_in_a_reg;
 reg [MULT_AMOUNT*CHANNEL_SIZE-1:0] mwf_to_mult_in_b_reg;
 reg [MULT_AMOUNT-1:0] mwf_to_mult_valid_reg;
 reg [MULT_AMOUNT-1:0] finished_mults;
 
 reg [31:0] mult_offset;
 reg [CHANNEL_SIZE-1:0] v1 [0:2*NUM_CHANNELS-1];
 reg [7:0] v1_index;
 
 wire [CHANNEL_SIZE-1:0] IMM [0:2*NUM_CHANNELS-1][0:2*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE-1:0] ryy_BETA_mult_result [0:2*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE-1:0] m1 [0:2*NUM_CHANNELS-1][0:2*NUM_CHANNELS-1];
 
 reg [CHANNEL_SIZE-1:0] w_mwf [0:2*NUM_CHANNELS-1];
 reg [CHANNEL_SIZE-1:0] coeff [0:NFFT-1][0:2*NUM_CHANNELS-1];
 
 reg update_word_reg;
 
 reg [CHANNEL_SIZE*2*NUM_CHANNELS-1:0] sample_coeff;
 reg coeff_sample_ack_reg;
 
 generate

        assign coeff_packed = sample_coeff;
        
    for (i=0; i<NUM_CHANNELS*NUM_CHANNELS; i=i+1) begin
        assign ry_w[i] = ry_packed_word[i*CHANNEL_SIZE +: CHANNEL_SIZE];
        assign rx_w[i] = rx_packed_word[i*CHANNEL_SIZE +: CHANNEL_SIZE];
    end
    for (i=0; i<NUM_CHANNELS; i=i+1) begin
        assign qLR[i] = qL_packed[i*(CHANNEL_SIZE/2) +: (CHANNEL_SIZE/2)];
        assign qLR[i+NUM_CHANNELS] = qL_packed[i*(CHANNEL_SIZE/2) +: (CHANNEL_SIZE/2)];
    end
    for (i=0; i<2*NUM_CHANNELS; i=i+1) begin
        for (j=0; j<2*NUM_CHANNELS; j=j+1) begin
            if (i == j) begin
                assign IMM[i][j] = {{(CHANNEL_SIZE/2){1'b0}}, 1'b1};
            end else begin 
                assign IMM[i][j] = 0;
            end
        end
    end
 endgenerate
 
 assign done_o = done_reg;
 assign ready_o = ready_reg;
 assign mwf_to_mult_valid = mwf_to_mult_valid_reg;
 assign mwf_to_mult_in_a = mwf_to_mult_in_a_reg;
 assign mwf_to_mult_in_b = mwf_to_mult_in_b_reg;
 assign sampe_addr_o = sample_counter;
 assign update_word_o = update_word_reg;
 assign coeff_sample_ack_o = coeff_sample_ack_reg;
 always @(posedge(clk)) begin
    done_reg <= 0;
    ready_reg <= 0;
    coeff_sample_ack_reg <= 0;
    if (rst) begin 
        update_word_reg <= 0;
        frame_valid_reg <= 0;
        sample_counter <= 0;
        for (channel_1=0; channel_1<NFFT; channel_1=channel_1+1) begin
            for (channel_2=0; channel_2<2*NUM_CHANNELS; channel_2=channel_2+1) begin
                if (channel_2 == 0 || channel_2 == 9) begin
                    coeff[channel_1][channel_2] <= 1;
                end else begin
                    coeff[channel_1][channel_2] <= 0;
                end
            end
        end
    end else begin
        case (state)
            READY: begin
                ready_reg <= 1;
                if (matrices_updated_i) begin
                    sample_counter <= 0;
                    update_word_reg <= 1;
                    state <= FETCH_M_WORDS;
                end else if (coeff_sample_request_i) begin
                    sample_counter <= coeff_sample_addr_i;
                    state <= SEND_COEFF_SAMPLE;
                end
            end
            SEND_COEFF_SAMPLE: begin
                for (channel_1=0; channel_1<2*NUM_CHANNELS; channel_1=channel_1+1) begin
                    sample_coeff[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= coeff[sample_counter][channel_1];
                end
                coeff_sample_ack_reg <= 1;
                ready_reg <= 1;
                state <= READY;
            end
            FETCH_M_WORDS: begin
                
                if (word_updated_i) begin 
                    update_word_reg <= 0;
                    state <= R_PREP;
                end
            end
            R_PREP: begin
                for (channel_1=0; channel_1<2*NUM_CHANNELS; channel_1=channel_1+1) begin
                    v1[channel_1] <= 0;
                    w_mwf[channel_1] <= 0;
                    for (channel_2=0; channel_2<2*NUM_CHANNELS; channel_2=channel_2+1) begin
                        if (channel_1 < NUM_CHANNELS) begin
                            if (channel_2 < NUM_CHANNELS) begin
                                rxx[channel_1][channel_2] <= rx_w[channel_1*NUM_CHANNELS*NUM_CHANNELS+channel_2];
                                ryy[channel_1][channel_2] <= ry_w[channel_1*NUM_CHANNELS*NUM_CHANNELS+channel_2];
                            end else begin 
                                rxx[channel_1][channel_2] <= 0;
                                ryy[channel_1][channel_2] <= 0;
                            end
                        end else begin 
                            if (channel_2 < NUM_CHANNELS) begin
                                rxx[channel_1][channel_2] <= 0;
                                ryy[channel_1][channel_2] <= 0;
                            end else begin 
                                rxx[channel_1][channel_2] <= rx_w[(channel_1-NUM_CHANNELS)*NUM_CHANNELS*NUM_CHANNELS+(channel_2-NUM_CHANNELS)];
                                ryy[channel_1][channel_2] <= ry_w[(channel_1-NUM_CHANNELS)*NUM_CHANNELS*NUM_CHANNELS+(channel_2-NUM_CHANNELS)];
                            end
                        end
                    end
                end
                state <= V1_MULT_1;
                finished_mults <= 0;
            end
            V1_MULT_1: begin
                for (channel_1=0; channel_1<2*NUM_CHANNELS; channel_1=channel_1+1) begin
                    mwf_to_mult_valid_reg[channel_1] <= 0;
                    if (mwf_mult_ready[channel_1]) begin
                        mwf_to_mult_in_a_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, BETA};
                        mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, qLR[channel_1]};
                        mwf_to_mult_valid_reg[channel_1] <= 1;
                    end
                    if (mwf_mult_done[channel_1]) begin
                        qLR_BETA_mult_result[channel_1] <= mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE];
                        
                        finished_mults[channel_1] <= 1;
                    end
                end
                if (&(finished_mults[2*NUM_CHANNELS-1:0])) begin
                    finished_mults <= 0;
                    mult_offset <= 0;
                    v1_index <= 0;
                    state <= V1_MULT_2;
                end
            end
            V1_MULT_2: begin
                for (channel_1=0; channel_1<REDUCED_MULT_AMOUNT; channel_1=channel_1+1) begin
                    mwf_to_mult_valid_reg[channel_1] <= 0;
                    if (mwf_mult_ready[channel_1]) begin
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            mwf_to_mult_in_a_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= qLR_BETA_mult_result[channel_1];
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= rxx[v1_index][channel_1];
                        end else begin
                            mwf_to_mult_in_a_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= qLR_BETA_mult_result[channel_1-2*NUM_CHANNELS];
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= rxx[v1_index+1][channel_1-2*NUM_CHANNELS];
                        end
                        mwf_to_mult_valid_reg[channel_1] <= 1;
                    end
                    if (mwf_mult_done[channel_1]) begin
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            v1[v1_index] <= v1[v1_index] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE];
                        end else begin
                            v1[v1_index+1] <= v1[v1_index+1] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE];
                        end
                        finished_mults[channel_1] <= 1;
                    end
                end
                if (&(finished_mults[REDUCED_MULT_AMOUNT-1:0])) begin
                    finished_mults <= 0;
                    v1_index <= v1_index + 2;
                    mult_offset <= mult_offset + REDUCED_MULT_AMOUNT;
                end
                if (mult_offset >= 4*NUM_CHANNELS*NUM_CHANNELS) begin
                    v1_index <= 0;
                    state <= M1_MULT; //talvez pular "V1_SUM" e ir direto pra M1_MULT
                end
            end
            V1_SUM: begin
                //provavelmente esse estado não é necessário
            end
            M1_MULT: begin
                for (channel_1=0; channel_1<REDUCED_MULT_AMOUNT; channel_1=channel_1+1) begin
                    mwf_to_mult_valid_reg[channel_1] <= 0;
                    if (mwf_mult_ready[channel_1]) begin
                        mwf_to_mult_valid_reg[channel_1] <= 1;
                        mwf_to_mult_in_a_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= {{(CHANNEL_SIZE/2){1'b0}}, BETA};
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= ryy[v1_index][channel_1];
                        end else begin
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= ryy[v1_index+1][channel_1-2*NUM_CHANNELS];
                        end
                    end
                    if (mwf_mult_done[channel_1]) begin
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            m1[v1_index][channel_1][CHANNEL_SIZE/2-1:0] <= IMM[v1_index][channel_1][CHANNEL_SIZE/2-1:0] - mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            m1[v1_index][channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= IMM[v1_index][channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] - mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end else begin
                            m1[v1_index+1][channel_1-2*NUM_CHANNELS][CHANNEL_SIZE/2-1:0] <= IMM[v1_index][channel_1-2*NUM_CHANNELS][CHANNEL_SIZE/2-1:0] - mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                            m1[v1_index+1][channel_1-2*NUM_CHANNELS][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= IMM[v1_index][channel_1-2*NUM_CHANNELS][CHANNEL_SIZE-1:CHANNEL_SIZE/2] - mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                        end
                        finished_mults[channel_1] <= 1;
                    end
                end
                if (&(finished_mults[REDUCED_MULT_AMOUNT-1:0])) begin
                    finished_mults <= 0;
                    v1_index <= v1_index + 2;
                    mult_offset <= mult_offset + REDUCED_MULT_AMOUNT;
                end
                if (mult_offset >= 4*NUM_CHANNELS*NUM_CHANNELS) begin
                    v1_index <= 0;
                    state <= MWF_CALC;
                end
            end
            MWF_CALC: begin
                for (channel_1=0; channel_1<REDUCED_MULT_AMOUNT; channel_1=channel_1+1) begin
                    mwf_to_mult_valid_reg[channel_1] <= 0;
                    if (mwf_mult_ready[channel_1]) begin
                        mwf_to_mult_valid_reg[channel_1] <= 1;
                        mwf_to_mult_in_a_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= coeff[sample_counter][channel_1];
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= m1[v1_index][channel_1];
                        end else begin
                            mwf_to_mult_in_b_reg[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE] <= m1[v1_index+1][channel_1-2*NUM_CHANNELS];
                        end
                    end
                    if (mwf_mult_done[channel_1]) begin
                        if (channel_1 < 2*NUM_CHANNELS) begin
                            if (channel_1 == 2*NUM_CHANNELS-1) begin
                                w_mwf[v1_index][CHANNEL_SIZE/2-1:0] <= w_mwf[v1_index][CHANNEL_SIZE/2-1:0] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2] + v1[v1_index][CHANNEL_SIZE/2-1:0];
                                w_mwf[v1_index][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= w_mwf[v1_index][CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + v1[v1_index][CHANNEL_SIZE-1:CHANNEL_SIZE/2];
                            end else begin 
                                w_mwf[v1_index][CHANNEL_SIZE/2-1:0] <= w_mwf[v1_index][CHANNEL_SIZE/2-1:0] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                                w_mwf[v1_index][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= w_mwf[v1_index][CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                            end
                        end else begin
                            if (channel_1 == 4*NUM_CHANNELS-1) begin
                                w_mwf[v1_index+1][CHANNEL_SIZE/2-1:0] <= w_mwf[v1_index+1][CHANNEL_SIZE/2-1:0] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2] + v1[v1_index+1][CHANNEL_SIZE/2-1:0];
                                w_mwf[v1_index+1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= w_mwf[v1_index+1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2] + v1[v1_index+1][CHANNEL_SIZE-1:CHANNEL_SIZE/2];
                            end else begin
                                w_mwf[v1_index+1][CHANNEL_SIZE/2-1:0] <= w_mwf[v1_index+1][CHANNEL_SIZE/2-1:0] + mwf_mult_result[channel_1*CHANNEL_SIZE +: CHANNEL_SIZE/2];
                                w_mwf[v1_index+1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= w_mwf[v1_index+1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] + mwf_mult_result[channel_1*CHANNEL_SIZE+CHANNEL_SIZE/2 +: CHANNEL_SIZE/2];
                            end
                            
                        end
                        finished_mults[channel_1] <= 1;
                    end
                end
                if (&(finished_mults[REDUCED_MULT_AMOUNT-1:0])) begin
                    finished_mults <= 0;
                    v1_index <= v1_index + 2;
                    mult_offset <= mult_offset + REDUCED_MULT_AMOUNT;
                end
                if (mult_offset >= 4*NUM_CHANNELS*NUM_CHANNELS) begin
                    v1_index <= 0;
                    state <= COEFF_UPDATE;
                end
            end
            COEFF_UPDATE: begin
                for (channel_1=0; channel_1<2*NUM_CHANNELS; channel_1=channel_1+1) begin
                    coeff[sample_counter][channel_1][CHANNEL_SIZE/2-1:0] <= w_mwf[channel_1][CHANNEL_SIZE/2-1:0];
                    coeff[sample_counter][channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= w_mwf[channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2];
                    if (sample_counter > 0 && sample_counter < (NFFT/2)) begin
                        coeff[NFFT-sample_counter][channel_1][CHANNEL_SIZE/2-1:0] <= w_mwf[channel_1][CHANNEL_SIZE/2-1:0];
                        coeff[NFFT-sample_counter][channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2] <= -$signed(w_mwf[channel_1][CHANNEL_SIZE-1:CHANNEL_SIZE/2]);
                    end
                    if (sample_counter == SAMPLE_USE-1) begin
                        ready_reg <= 1;
                        state <= READY;
                        sample_counter <= 0;
                        done_reg <= 1;
                    end else begin 
                        sample_counter <= sample_counter + 1;
                        update_word_reg <= 1;
                        state <= FETCH_M_WORDS;
                    end
                end
            end
            default: begin
                ready_reg <= 1;
                state <= READY;
                end
        endcase
    end
 end
endmodule
