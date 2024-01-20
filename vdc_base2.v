module vdc_base2 #(parameter width=13) (clk, reset, vdcout);
input clk,reset;
output [width-1:0] vdcout;
reg [width-1:0] outcount;
always@(posedge clk)
begin
    if(reset)
        outcount <= 0;
    else
        outcount <= outcount + 1;
end
assign vdcout [width-1:0] = {outcount[0],outcount[1],outcount[2],outcount[3],outcount[4],outcount[5],outcount[6],outcount[7],outcount[8],outcount[9],outcount[10],outcount[11],outcount[12]};
endmodule
