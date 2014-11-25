#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <sys/time.h>

#include <mpi.h> //mpi
#define _REENTRANT //mpi

#define C 2.0
#define R 2.0
//#define I 5.0
#define T2C 10.0 //capasity - time relation
#define CUR 2.0 // currency in one of nodes
#define CUR_NOD 0// number of node with curr!=0

double *volt_cur;//mass inside time step
double *volt_last;//mass inside last time step
double *volt_last_hard;//mass inside time step
double *flow;//mass of currencies
pthread_t* part_by_thread;
int threads_num;
int nodes_len;
int nodes_wid;
int nodes_num;
int time_num;
int nodes_per_thread;
double ht;
pthread_barrier_t bp;

//timing
struct timeval tv1,tv2,dtv;
struct timezone tz;

void time_start() { gettimeofday(&tv1, &tz); }
long time_stop()
{ gettimeofday(&tv2, &tz);
    dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
    dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
    if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
    return dtv.tv_sec*1000+dtv.tv_usec/1000;
}

void visualiser(void)
{
    int i,j;
    FILE *fp;
    fp = fopen ("gnulab.txt","a");
    fprintf(fp, "plot '-' using 1:2 smooth csplines\n");
    //for(j=0;j<2;j++)
    for(i=0; i<nodes_num; i++)
    {	fprintf(fp, "%d ", i);
        fprintf(fp, "%f\n", volt_last[i]);
        if(!((i+1)%nodes_wid))
        {
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "e\npause 0.01\n\n" );
    fclose(fp);
}

int getFirstById(int id)
{	
    return id*nodes_per_thread;
}

int getLastById(pthread_t id)
{
    if(id==threads_num-1)
        return nodes_num;
    return id*nodes_per_thread + nodes_per_thread;
}

char whatIs(int i)
{

    if(!i)
    {
        //leftside
        return 1;
    }
    else if (i==nodes_num-1)
    {
        //rightside
        return 2;
    }
    else
        //inside
        return 3;

}

void compute(int first, int last, int myrank)
{
    int t,i,u,timer=0;
    MPI_Status stat;
    double * volt_cur_mini = malloc(sizeof(double));
    volt_cur_mini = (double *)calloc(last-first, sizeof(double));
    for(t=0; t < time_num; t++)
    {	
        //if(t==0)time_start();//////////
        for(i=first; i<last; i++)
        {	
            char what = whatIs(i);
            switch(what)
            {
                case 1: //leftside
                    volt_cur[i] = volt_last[i] +
                        (ht/C)*(flow[i] + (1.0/R)*
                                (-2*volt_last[i] + volt_last[i+1]));
                    break;
                case 2: //rightside
                    volt_cur[i] = volt_last[i] +
                        (ht/C)*(flow[i] + (1.0/R)*
                                (-2*volt_last[i] + volt_last[i-1]));
                    break;
                case 3: //inside
                    volt_cur[i] = volt_last[i] +
                        (ht/C)*(flow[i] + (1.0/R)*
                                (-2.0*volt_last[i] + volt_last[i+1]+ 
                                 volt_last[i-1]));
            }

        }
        //barrier !!!!!
        MPI_Barrier(MPI_COMM_WORLD);


        //swapping

        //for(i = 0; i<(last-first); volt_cur_mini[i] = volt_cur[i+first],i++);
        for(i = first; i<(last); volt_last[i] = volt_cur[i],i++) {
        }
    


        MPI_Barrier(MPI_COMM_WORLD);

        if(threads_num!=1)
        {
            if(myrank==threads_num-1)
            {//rightside
                MPI_Send((void*)&(volt_last[first]),1,
                        MPI_DOUBLE,myrank-1,1,MPI_COMM_WORLD);
                MPI_Recv((void*)&(volt_last[first-1]),1,
                        MPI_DOUBLE,myrank-1,MPI_ANY_TAG,
                        MPI_COMM_WORLD,&stat);
            }
            else if(!myrank)
            {//leftside
                MPI_Send((void*)&(volt_last[last-1]),1,
                        MPI_DOUBLE,myrank+1,1,MPI_COMM_WORLD);
                MPI_Recv((void*)&(volt_last[last]),1,
                        MPI_DOUBLE,myrank+1,MPI_ANY_TAG,
                        MPI_COMM_WORLD,&stat);
            }
            else
            {
                MPI_Send((void*)&(volt_last[last-1]),1,
                        MPI_DOUBLE,myrank+1,1,MPI_COMM_WORLD);
                MPI_Send((void*)&(volt_last[first]),1,
                        MPI_DOUBLE,myrank-1,1,MPI_COMM_WORLD);
                MPI_Recv((void*)&(volt_last[last]),1,
                        MPI_DOUBLE,myrank+1,MPI_ANY_TAG,
                        MPI_COMM_WORLD,&stat);
                //volt_last[last]=*volt_cur_mini;
                MPI_Recv((void*)&(volt_last[first-1]),1,
                        MPI_DOUBLE,myrank-1,MPI_ANY_TAG,
                        MPI_COMM_WORLD,&stat);
                //volt_last[first-1]=*volt_cur_mini;
            }
        }

        // MPI_Bcast((void *)volt_last, nodes_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //visualise volt_last
        MPI_Barrier(MPI_COMM_WORLD);
        //barrier !!!!!
        //if(t==0)fprintf(stderr,"Time: %ld, myrank is %d\n", time_stop(),myrank);
        MPI_Gather((void *)&(volt_last[first]), last-first, MPI_DOUBLE,
                (void *)&(volt_last[0]), last-first, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!myrank) {visualiser();}
    }

}


int main(int argc, char const *argv[])
{
    nodes_len = 1;
    nodes_wid = atoi(argv[1]);
    //threads_num = atoi(argv[2]);
    int m_time = atoi(argv[2]);
    ht = (C/T2C) * 1e-1;

    int myrank, total, i;
    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &total);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    printf ("Total=%d, rank=%d\n", total, myrank);

    threads_num = total;
    nodes_num = nodes_wid*nodes_len;
    if(nodes_num%threads_num)
        nodes_per_thread = (int)((double)nodes_num/(double)threads_num)+1;
    else
        nodes_per_thread = (int)((double)nodes_num/(double)threads_num);
    time_num = (int)(m_time/ht);
    
    volt_cur = (double *)calloc(nodes_num, sizeof(double));
    for(i = 0; i<nodes_num; volt_cur[i++] = 0.0);
    volt_last = (double *)calloc(nodes_num, sizeof(double));
    for(i = 0; i<nodes_num; volt_last[i++] = 0.0);
    //volt_last[0]=100;
    //volt_last[nodes_num-1]=50;

    //allocate currency memory
    flow = (double *)calloc(nodes_wid, sizeof(double));
    for(i = 0; i < nodes_wid; flow[i++] = 0.0);
    flow[nodes_wid / 2]=CUR;
    if(!myrank)
    {
        //super hard array
        volt_last_hard = (double *)calloc(nodes_num, sizeof(double));
        for(i = 0; i<nodes_num; volt_last_hard[i++] = 0.0);
        //time steps number
        if(threads_num>nodes_num)
        {
            perror("threads num can't be grater then nodes num");
            exit(2);
        }

        //allocate voltage memory

        //prepare file for gnuplt
        FILE *fp;
        fp = fopen ("gnulab.txt","w");
        fprintf(fp, "set cbrange [0:20]\n");
        fprintf(fp, "set cblabel \"Score\"\n");
        fprintf(fp, "set xrange [0:%d]\n", nodes_wid-1);
        fprintf(fp, "set yrange [0:%d]\n", 50);
        fprintf(fp, "set view map\n");
        //fprintf(fp, "set size ratio 1\n"); ?????????????????????????
        fclose(fp);

        visualiser();	
        //pthread_barrier_init(&bp, NULL, threads_num);
        //create threads and pass they to function 
        time_start();
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    compute(getFirstById(myrank),getLastById(myrank),myrank);

    if(!myrank)
    {	
        fprintf(stderr,"Time: %ld\n", time_stop());
        system("gnuplot -persist \"gnulab.txt\""); //Ñ‡Ñ‚Ð¾Ð±Ñ‹ Ð¿Ñ€Ð¾Ð²ÐµÑ€Ð¸Ñ‚ÑŒ Ñ€Ð°Ð±Ð¾Ñ‚Ñƒ
        
    }

    // pthread_barrier_destroy(&bp);
    MPI_Finalize();
    exit(0);
}

