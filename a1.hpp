/*  YOUR_FIRST_NAME Parush 
 *  YOUR_LAST_NAME Garg 
 *  YOUR_UBIT_NAME parushga
 */

#ifndef A1_HPP
#define A1_HPP
#define barrier() MPI_Barrier(comm)

#include <unistd.h>
#include <vector>
#include <mpi.h>

int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    // ...
    int len, com_col, com_row, rank, max;
    MPI_Comm_rank(comm,&rank);

    len = n/q;
    int col = rank % q;
    int row = rank / q;

    std::vector<int> p; // Initial P Vector
    std::vector<int> p_prime(len*len,0); // P vector after opportunistic pointer jumping
    std::vector<int> pTreeHanging; // P vector after Tree hanging
    std::vector<int> m(len*len,0); // Helper matrix
    std::vector<int> c(len*len,0); // Helper matrix
    std::vector<int> vectorQ(len*len,0); // Helper matrix
    std::vector<int> fullPVector(n*n,0); 
    //int diff = 1; // if there is difference in initial p vector and p tree hanging vector. Initialised as true

    // j*len + i -> row wise
    // i*len + j -> column wise

    // Initialise P vector if A[i] contains element 1
    for(int i=0; i< static_cast<int>(A.size()); i += len)  
    {
        for (int j = i; j< len + i; j++)
        {
            if(static_cast<int>(A[j]) == 1)
            {
                int k = (j/len) + (row*len);
                p.push_back(k);
            }
            else 
            {
                p.push_back(0);
            }
        }
    }

    //To get sync with all ranks
    barrier();
    sleep(rank);

    // while(diff)
    // {
        diff = 0;
        std::cout<<"P :\n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << p[(i * len) + j] << " ";
            std::cout << std::endl;
        }

        barrier();
        //Finding max sequentially
        for (int i = 0; i < len; ++i) 
        {
            max=0;
                for (int j = 0; j < len; ++j)
                {
                    if(max<p[(j*len)+i])
                        max = p[(j*len)+i];
                }
                for (int j = 0; j < len; ++j)
                {
                        p[(j*len)+i] = max;
                }
        }

        barrier();

        //Finding vector P by doing Allreduce to M and copying M to P - Colwise
        MPI_Comm_split(comm,col,rank,&com_col);
        MPI_Allreduce(p.data(),m.data(),(len*len),MPI_INT,MPI_MAX,com_col);
        MPI_Comm_free(&com_col);

        barrier();

        for(int i=0;i<static_cast<int>(m.size());i++)
            p[i] = m[i];

        barrier();
        sleep(rank);

        std::cout<<"PInit : \n";                // PInit which is first step in Algorithm
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << p[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        //Making helper matrix M to compute vectorQ 
        for(int i=0;i<static_cast<int>(A.size());++i)
        {
            if(static_cast<int>(A[i])==1)
                m[i]=p[i];
            else
                m[i] = 0;
        }

        barrier();
        sleep(rank);

        std::cout<<"M : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << m[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        //Finding max sequentially
        for (int i = 0; i < len; ++i) {
            max=0;
            for (int j = 0; j < len; ++j)
            {
                if(max < m[i*len+j])
                    max = m[i*len+j];
            }
            for (int j = 0; j < len; ++j)
            {
                m[i*len+j] = max;
            }
        }

        barrier();

        //Computing vector Q - Row wise Allreduce
        MPI_Comm_split(comm,row,rank,&com_row);
        MPI_Allreduce(m.data(),vectorQ.data(),(len*len),MPI_INT,MPI_MAX,com_row);
        MPI_Comm_free(&com_row);

        barrier();
        sleep(rank);

        std::cout<<"vectorQ : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << vectorQ[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        //Computing vector M now to find the opportunistic pointer jumping vector P
        for(int i=0;i<len;i++)
        {
            for(int j=0;j<len;j++)
            {
            if(vectorQ[i*len+j]==j+(len*col))
            {
                m[i*len+j] = p[i*len+j];
            }
            else
            {
                m[i*len+j] = 0;
            }
            }
        }

        /*
        for(int i = 0; i < len ; i++)
        {
            for(int j = 0; j < len; j++)
            {
                if(c[j+(len * i)] == j + ((rank / q) * len))
                    m[j+(len * i)] = p[j+(len * i)];
                else 
                    m[j+(len * i)] = 0;
                // if((i + ((rank / q) * len)) == 6)
                //     std::cout<<"P: "<<p[j+(len * i)]<<"::"<<"PPrime: "<<p_prime[j+(len * i)]<<std::endl;
            }
        } */



        barrier();
        sleep(rank);

        std::cout<<"M1 : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << m[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        //Finding max sequentially
        for (int i = 0; i < len; ++i) {
            max=0;
            for (int j = 0; j < len; ++j)
            {
                if(max<m[i*len+j])
                    max = m[i*len+j];
            }
            for (int j = 0; j < len; ++j)
            {
                m[i*len+j] = max;
            }

        }
        barrier();

        //Finding P_Prime which is my vector after opportunistic pointer jumping
        MPI_Comm_split(comm,row,rank,&com_row);
        MPI_Allreduce(m.data(),p_prime.data(),(len*len),MPI_INT,MPI_MAX,com_row);
        MPI_Comm_free(&com_row);

        barrier();
        sleep(rank);

        std::cout<<"PPrime : \n";                               // PPrime which is second step in Algorithm
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << p_prime[i * len + j] << " ";
            std::cout << std::endl;
        }

        // Tree Hanging step

        sleep(rank);
        barrier();

        for(int i = 0; i < len ; i++)
        {
            for(int j = 0; j < len; j++)
            {
                if(p[j+(len * i)] == i + ((rank / q) * len))
                    m[j+(len * i)] = p_prime[j+(len * i)];
                else 
                    m[j+(len * i)] = 0;
                // if((i + ((rank / q) * len)) == 6)
                //     std::cout<<"P: "<<p[j+(len * i)]<<"::"<<"PPrime: "<<p_prime[j+(len * i)]<<std::endl;
            }
        }

        barrier();
        sleep(rank);

        std::cout<<"MTreeHanging : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << m[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();
        sleep(rank);

        //Finding max sequentially
        for (int i = 0; i < len; ++i) {
            max=0;
            for (int j = 0; j < len; ++j)
            {
                if(max<m[i*len+j])
                    max = m[i*len+j];
            }
            for (int j = 0; j < len; ++j)
            {
                m[i*len+j] = max;
            }
        }

        barrier();

        //Finding C vector which is Q in the example explained in the class
        MPI_Comm_split(comm,row,rank,&com_row);
        MPI_Allreduce(m.data(),c.data(),(len*len),MPI_INT,MPI_MAX,com_row);
        MPI_Comm_free(&com_row);

        barrier();
        sleep(rank);

        std::cout<<"CTreeHanging : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << c[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        //Finding pTreeHanging vector after Tree Hanging step.
        for(int i=0 ; i< len; i++)
        {
            for(int j = 0; j< len; j++)
            {
            if(p_prime[i*len + j] < c[i*len + j]) // maximum between c[] and p_prime[]
                pTreeHanging.push_back(c[i*len + j]);
            else 
                pTreeHanging.push_back(p_prime[i*len + j]);
            }
        }

        barrier();
        sleep(rank);

        std::cout<<"PTreeHanging : \n";                     // PTreeHanging which is third step in Algorithm
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << pTreeHanging[i * len + j] << " ";
            std::cout << std::endl;
        }

        barrier();

        // for (int i = 0; i < len; ++i) {
        //     for (int j = 0; j < len; ++j)
        //     std::cout << pTreeHanging[i * len + j] << " ";
        //     std::cout << std::endl;
        // }
        
        
        /*
        //MPI_Allreduce(m.data(),p_prime.data(),(len*len),MPI_INT,MPI_MAX,com_row);

        //MPI_Gather(pTreeHanging.data(), n*n, MPI_INT, fullPVector.data(), n*n, MPI_INT, 0 ,MPI_Comm comm);

        //Checking if PPrime after round1 is same as PInitialised
        //std::cout<<"Rank :"<<rank<<"\n";
        for(int i=0 ; i< len; i++)
        {
            for(int j = 0; j< len; j++)
            {

                //std::cout<<"PVector : "<<p[i*len + j]<<"pTreeHangingVector : "<<pTreeHanging[j*len + i]<<std::endl;
                if(p[i*len + j] != pTreeHanging[j*len + i])
                {
                    //std::cout<<"M1 : "<<p[i*len + j]<<"M2 : "<<pTreeHanging[j*len + i]<<std::endl;
                    diff = 1;
                    break;
                }
                else 
                {
                    diff = 0;
                }
            }
        }
        barrier();
        //sleep(rank);
        int globaldiff;
        //MPI_Gather(diff, 1, MPI_BOOL, fullPVector.data(), n*n, MPI_INT, 0 ,MPI_Comm comm);
        MPI_Comm_split(comm,row,rank,&com_row);
        //MPI_Allreduce(diff, globaldiff, 1, MPI_INT, MPI_SUM, com_row);
        MPI_Gather(p.data(), len*len, MPI_INT, fullPVector.data(), len*len, MPI_INT, 0 ,com_row);
        MPI_Comm_free(&com_row);
        barrier();
        sleep(rank);
        std::cout<<"Final P Vector : \n";
        std::cout<<"Rank :"<<rank<<"\n";
        for (int i = 0; i < len; ++i) {
            for (int j = 0; j < len; ++j)
            std::cout << fullPVector[i * len + j] << " ";
            std::cout << std::endl;
        }
        //std::cout<<"Are Vector P's different?  "<<diff<<std::endl;
        p = pTreeHanging;
        */
    //}
    return -1;
} // connected_components

#endif // A1_HPP

//module load intel-mpi/2017.0.1
