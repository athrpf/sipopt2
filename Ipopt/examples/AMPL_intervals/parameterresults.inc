
#printf "%5.5f\n",xAL_init>xAL.dat;
printf "" > xAL.dat;
printf {i in int, j in fe,k in cp}:
"%5.5f\n",xAL[i,j,k]>>xAL.dat;

#printf "%5.5f\n",xAU_init>xAU.dat;
printf "" > xAU.dat;
printf {i in int, j in fe,k in cp}:
"%5.5f\n",xAU[i,j,k]>>xAU.dat;

#printf "%5.5f\n",xBL_init>xBL.dat;
printf "" > xBL.dat;
printf {i in int, j in fe,k in cp}:
"%5.5f\n",xBL[i,j,k]>>xBL.dat;

#printf "%5.5f\n",xBU_init>xBU.dat;
printf "" > xBU.dat;
printf {i in int, j in fe,k in cp}:
"%5.5f\n",xBU[i,j,k]>>xBU.dat;

printf "%5.5f\n", time > time.dat;
