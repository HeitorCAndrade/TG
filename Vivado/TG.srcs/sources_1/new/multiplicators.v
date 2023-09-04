`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 18.08.2023 10:30:45
// Design Name: 
// Module Name: multiplicators
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


module multiplicators #(
    parameter DECIMAL_BITS = 20,
    parameter CHANNEL_SIZE = 64,
    parameter MULT_AMOUNT = 32
)
(
    input wire clk,
    input wire rst,
    input wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] input_a,
    input wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] input_b,
    input wire [MULT_AMOUNT-1:0] valid,
    input wire [MULT_AMOUNT-1:0] real_mult,
    output wire [MULT_AMOUNT*CHANNEL_SIZE-1:0] result,
    output wire [MULT_AMOUNT-1:0] ready,
    output wire [MULT_AMOUNT-1:0] done,
    output wire [MULT_AMOUNT-1:0] ovflw
);

localparam READY = 0,
           MULT = 1;

genvar i;
integer index;

reg state;

wire [CHANNEL_SIZE-1:0] input_a_w [0:MULT_AMOUNT-1];
wire [CHANNEL_SIZE-1:0] input_b_w [0:MULT_AMOUNT-1];
wire [CHANNEL_SIZE-1:0] cmp_result [0:MULT_AMOUNT-1];
//reg cmp_valid [0:MULT_AMOUNT-1];
wire cmp_ready [0:MULT_AMOUNT-1];
wire cmp_done [0:MULT_AMOUNT-1];
wire cmp_ovfl [0:MULT_AMOUNT-1];

generate
    for (i=0; i<MULT_AMOUNT; i=i+1) begin
        assign result[i*CHANNEL_SIZE +: CHANNEL_SIZE] = cmp_result[i];
        assign ovflw[i] = cmp_ovfl[i];
        assign ready[i] = cmp_ready[i];
        assign done[i] = cmp_done[i];
        assign input_a_w[i] = input_a[i*CHANNEL_SIZE +: CHANNEL_SIZE];
        assign input_b_w[i] = input_b[i*CHANNEL_SIZE +: CHANNEL_SIZE];
        complex_mult #(
            .NUMBER_WIDTH(CHANNEL_SIZE/2),
            .DECIMAL_BITS(DECIMAL_BITS)
        ) complex_mult (
            .clk(clk),
            .rst(rst),
            .input_a(input_a_w[i]),
            .input_b(input_b_w[i]),
            .in_valid(valid[i]),
            .result(cmp_result[i]),
            .ready(cmp_ready[i]),
            .done(cmp_done[i]),
            .ovfl_flag(cmp_ovfl[i])
        ); 
    end
endgenerate
endmodule
