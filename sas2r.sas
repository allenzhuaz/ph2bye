/* save SAS dataset in trasport format */
libname out xport 'c:/mydata.xpt';
  data out.mydata;
  set sasuser.mydata;
run;
