// Down Counter by Mehran.
module down_counter #(parameter WIDTH=8) (clk, data_in, zero);

input clk;
//input rst;
input [7:0] data_in;
//output reg[WIDTH-1:0] out_count;
output wire zero;

reg [7:0] out_count;

//reg [WIDTH-1:0] counter_down;

// down counter
//assign out_count = data_in;
//assign zero =0;
always @(posedge clk)
begin
	if (zero)
		out_count <= data_in;
	else
		//out_count <= data_in;
		//zero <= 0;
	
		out_count <= out_count - 1;
	//end
 
	 //if (out_count == 0)
		//zero <= 1;
	//end 
end

//assign zero = ~out_count[0]&~out_count[1]&~out_count[2]&~out_count[3]&~out_count[4]&~out_count[5]&~out_count[6]&~out_count[7];
//assign zero = ~(out_count[0]|out_count[1]|out_count[2]|out_count[3]|out_count[4]|out_count[5]|out_count[6]|out_count[7]);
assign zero = ~|(out_count);
endmodule