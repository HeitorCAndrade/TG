`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 05.08.2023 11:29:30
// Design Name: 
// Module Name: complex_mult_tb
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


module complex_mult_tb;
    localparam NUM_TESTS = 4;
    localparam NUM_TESTS_LOG = 2;
    localparam NUMBER_WIDTH = 32;
    localparam DECIMAL_BITS = 20;
    
    integer num_accepted = 0;
    integer num_acc_overflow = 0;
    integer num_tests = 0;
    
    reg clk;
    reg rst;
    
    reg [2*2*NUMBER_WIDTH-1:0] mem_in [0:NUM_TESTS-1];
    reg [2*NUMBER_WIDTH:0] mem_out [0:NUM_TESTS-1];
    reg [NUM_TESTS_LOG-1:0] addr_in = 0;
    reg [NUM_TESTS_LOG-1:0] addr_out = 0;
    
    reg [2*NUMBER_WIDTH-1:0] input_a_reg = 0;
    reg [2*NUMBER_WIDTH-1:0] input_b_reg = 0;
    reg input_valid = 0;
    wire overflow_w;
    wire ready_w;
    wire done_w;
    wire [2*NUMBER_WIDTH-1:0] output_w;
    reg [2*NUMBER_WIDTH-1:0] output_reg = 0;
    reg [2*NUMBER_WIDTH-1:0] expected_out;
    reg expected_overflow;
    reg compare_out = 0;
    
    complex_mult # (
        .NUMBER_WIDTH(NUMBER_WIDTH),
        .DECIMAL_BITS(DECIMAL_BITS)
    ) dut(
        .clk(clk),
        .rst(rst),
        .input_a(input_a_reg),
        .input_b(input_b_reg),
        .in_valid(input_valid),
        .result(output_w),
        .ready(ready_w),
        .done(done_w),
        .ovfl_flag(overflow_w)
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
    
    initial begin
        $readmemb("tb_complex_mult_in.txt", mem_in); 
        $readmemb("tb_complex_mult_out.txt", mem_out); 
        $display("###################### LOADED MEM #######################");
    end
    
    always @(posedge(clk)) begin 
        //output_reg <= 0;
        compare_out <= 0;
        input_valid <= 0;
        input_a_reg <= mem_in[addr_in][0 +: 2*NUMBER_WIDTH];
        input_b_reg <= mem_in[addr_in][2*NUMBER_WIDTH +: 2*NUMBER_WIDTH];
        expected_out <= mem_out[addr_out][0 +: 2*NUMBER_WIDTH];
        expected_overflow <= mem_out[addr_out][2*NUMBER_WIDTH];
        if (ready_w) begin 
            input_valid <= 1;
            addr_in <= addr_in + 1;
        end
        if (done_w) begin
            output_reg <= output_w;
            compare_out <= 1;
            addr_out <= addr_out + 1;
        end
    end
    
    always @(compare_out) begin 
        if (compare_out) begin 
            num_tests = num_tests + 1;
            if (expected_out == output_reg) begin 
                num_accepted = num_accepted + 1;
                $display("test number %d accepted!", num_tests);
            end else begin
                $display("test number %d failed!", num_tests);
            end
            if (expected_overflow == overflow_w) begin 
                num_acc_overflow = num_acc_overflow + 1;
                $display("expected overflow");
            end else begin 
                if (expected_overflow) begin 
                    $display("oveflow on golden model");
                end else begin 
                    $display("oveflow on dut");
                end
            end
            if (num_tests == NUM_TESTS) begin 
                $display("##################TEST ENDED####################");
                if (num_tests == num_accepted) begin 
                    $display("TEST PASSED! ALL INSTANCES (%d) CORRECT", num_accepted);
                end else begin 
                    $display("TEST FAILED! %d OUT OF %d WRONG", num_tests-num_accepted, num_tests);
                end
                $finish;
            end
        end
    end

endmodule
