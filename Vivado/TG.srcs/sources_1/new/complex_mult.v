`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 04.08.2023 15:45:38
// Design Name: 
// Module Name: complex_mult
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


module complex_mult #(
    parameter NUMBER_WIDTH = 32,
    parameter DECIMAL_BITS = 20
)
(
    input wire clk,
    input wire rst,
    input wire [2*NUMBER_WIDTH-1:0] input_a,
    input wire [2*NUMBER_WIDTH-1:0] input_b,
    input wire in_valid,
    output wire [2*NUMBER_WIDTH-1:0] result,
    output wire ready,
    output wire done,
    output wire ovfl_flag
    );
    
    localparam READY = 0,
               MULT = 1;
    
    reg signed [NUMBER_WIDTH-1:0] MAX_VAL = {1'h0, {NUMBER_WIDTH-1{1'h1}}};
    reg signed [NUMBER_WIDTH-1:0] MIN_VAL = {1'h1, {NUMBER_WIDTH-1{1'h0}}};
    
    reg signed [2*NUMBER_WIDTH:0] MAX_VAL_CMP = 65'h00007FFFFFFF00004;
    reg signed [2*NUMBER_WIDTH:0] MIN_VAL_CMP = 65'h1FFF8000000000000;
    
    reg state;
    reg ready_reg;
    reg done_reg;
    
    reg signed [NUMBER_WIDTH-1:0] input_a_re;
    reg signed [NUMBER_WIDTH-1:0] input_a_im;
    reg signed [NUMBER_WIDTH-1:0] input_b_re;
    reg signed [NUMBER_WIDTH-1:0] input_b_im;
    reg signed [2*NUMBER_WIDTH-1:0] re_mult[0:1];
    reg [2*NUMBER_WIDTH:0] re_add;
    reg signed [2*NUMBER_WIDTH-1:0] imag_mult[0:1];
    reg [2*NUMBER_WIDTH:0] imag_add;
    reg [NUMBER_WIDTH-1:0] result_re;
    reg [NUMBER_WIDTH-1:0] result_imag;
    reg overflow_re;
    reg overflow_imag;
    reg [1:0] pipe_stage;
    
    wire expected_re;
    wire expected_imag;
    wire sign_re;
    wire sign_imag;
    
    assign ready = ready_reg;
    assign done = done_reg;
    
    assign expected_re = re_add[2*NUMBER_WIDTH];
    assign expected_imag = imag_add[2*NUMBER_WIDTH];
    
    assign sign_re = re_add[DECIMAL_BITS+NUMBER_WIDTH-1];
    assign sign_imag = imag_add[DECIMAL_BITS+NUMBER_WIDTH-1];
    
    assign ovfl_flag = overflow_re | overflow_imag;
    assign result = {result_imag, result_re};
    
    always @* begin 
        re_add = re_mult[0] - re_mult[1];
        imag_add = imag_mult[0] + imag_mult[1];
    end
    
    always @(posedge(clk)) begin 
        done_reg <= 0;
        if (rst) begin 
            state <= READY;
            ready_reg <= 0;
            input_a_re <= 0;
            input_a_im <= 0;
            input_b_re <= 0;
            input_b_im <= 0;
            pipe_stage <= 0;
            result_imag <= 0;
            result_re <= 0;
        end else begin 
            case(state)
                READY: begin 
                    ready_reg <= 1;
                    if (in_valid) begin 
                        input_a_re <= input_a[NUMBER_WIDTH-1:0];
                        input_a_im <= input_a[2*NUMBER_WIDTH-1:NUMBER_WIDTH];
                        input_b_re <= input_b[NUMBER_WIDTH-1:0];
                        input_b_im <= input_b[2*NUMBER_WIDTH-1:NUMBER_WIDTH];
                        state <= MULT;
                        ready_reg <= 0;
                        pipe_stage <= 1;
                    end
                end 
                MULT: begin 
                    pipe_stage <= pipe_stage << 1;
                    re_mult[0] <= input_a_re * input_b_re; //real * real
                    re_mult[1] <= input_a_im * input_b_im; //imag * imag
                
                    imag_mult[0] <= input_a_re * input_b_im; //real * imag
                    imag_mult[1] <= input_a_im * input_b_re; //imag * real
                    
                    if ((re_add < MIN_VAL_CMP) && re_add[2*NUMBER_WIDTH]) begin 
                        result_re <= MIN_VAL;
                        overflow_re <= 1;
                    end else begin
                        if ((re_add > MAX_VAL_CMP) && ~re_add[2*NUMBER_WIDTH]) begin 
                            result_re <= MAX_VAL;
                            overflow_re <= 1;
                        end else begin 
                            result_re <=  re_add[DECIMAL_BITS +: NUMBER_WIDTH];
                            overflow_re <= 0;
                        end  
                    end
                    if ((imag_add < MIN_VAL_CMP) && imag_add[2*NUMBER_WIDTH]) begin 
                        result_imag <= MIN_VAL;
                        overflow_imag <= 1;
                    end else begin
                        if ((imag_add > MAX_VAL_CMP) && ~imag_add[2*NUMBER_WIDTH]) begin 
                            result_imag <= MAX_VAL;
                            overflow_imag <= 1;
                        end else begin 
                            result_imag <=  imag_add[DECIMAL_BITS +: NUMBER_WIDTH];
                            overflow_imag <= 0;
                        end  
                    end
                    if (pipe_stage[1] == 1) begin 
                        done_reg <= 1;
                        ready_reg <= 1;
                        state <= READY;
                    end
                end 
            endcase
        end
    end
    
endmodule
