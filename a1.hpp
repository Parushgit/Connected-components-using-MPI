/*  YOUR_FIRST_NAME
 *  YOUR_LAST_NAME
 *  YOUR_UBIT_NAME
 */

#ifndef A1_HPP
#define A1_HPP

#include <unistd.h>
#include <vector>
#include <mpi.h>
int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    // ...
    int length,rank,size,max,c1;
    MPI_Comm_rank(comm,&rank);
    int col = rank % q;
    int row = rank / q;
    length = n/q;
    int length1 = n/q;
    int length2 = n/q;
    std::vector<int> p;
    std::vector<int> p_prime(length2*length2,0);
    std::vector<int> v;
    std::vector<int> pTreeHanging;
    std::vector<int> temp;
    std::vector<int> temp1(length,0);
    std::vector<int> m(length*length,0);
    std::vector<int> c(length1*length1,0);
    for(int i=0;i<A.size();++i)
    {
      if(static_cast<int>(A[i]) == 1)
      {
        p.push_back((i/length)+(row*length));
      }
      else
      {
        p.push_back(0);
      }
    }
    MPI_Barrier(comm);
    sleep(rank);
    std::cout<<"P :\n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << p[i * length + j] << " ";
        std::cout << std::endl;
    }
    MPI_Barrier(comm);
    for (int i = 0; i < length; ++i) {
        max=0;
        for (int j = 0; j < length; ++j)

          if(max<p[j*length+i])
            max = p[j*length+i];
        for (int j = 0; j < length; ++j)
        {
            p[j*length+i] = max;
        }

    }
    MPI_Barrier(comm);
    MPI_Comm_split(comm,col,rank,&c1);
    MPI_Allreduce(p.data(),m.data(),(length*length),MPI_INT,MPI_MAX,c1);
    MPI_Comm_free(&c1);
    MPI_Barrier(comm);
    for(int i=0;i<m.size();i++)
    {
      p[i] = m[i];
    }
    MPI_Barrier(comm);
    sleep(rank);
    std::cout<<"P1 : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << p[i * length + j] << " ";
        std::cout << std::endl;
    }
    MPI_Barrier(comm);
    for(int i=0;i<A.size();++i)
    {
      if(static_cast<int>(A[i])==1)
      {
        m[i]=p[i];
      }
      else
      {
        m[i] = 0;
      }
    }
    MPI_Barrier(comm);
    sleep(rank);
    std::cout<<"M : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << m[i * length + j] << " ";
        std::cout << std::endl;
    }
    MPI_Barrier(comm);
    for (int i = 0; i < length; ++i) {
        max=0;
        for (int j = 0; j < length; ++j)
          if(max<m[i*length+j])
            max = m[i*length+j];
        for (int j = 0; j < length; ++j)
        {
            m[i*length+j] = max;
        }

    }
    MPI_Barrier(comm);
    MPI_Comm_split(comm,row,rank,&c1);
    MPI_Allreduce(m.data(),c.data(),(length*length),MPI_INT,MPI_MAX,c1);
    MPI_Barrier(comm);
    sleep(rank);
    std::cout<<"C : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << c[i * length + j] << " ";
        std::cout << std::endl;
    }
    MPI_Barrier(comm);
    for(int i=0;i<length;i++)
    {
      for(int j=0;j<length;j++)
      {
        if(c[i*length+j]==j+(length*col))
        {
          m[i*length+j] = p[i*length+j];
        }
        else
        {
          m[i*length+j] = 0;
        }
      }
    }
    MPI_Barrier(comm);
    sleep(rank);
    std::cout<<"M1 : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << m[i * length + j] << " ";
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
    for (int i = 0; i < length; ++i) {
        max=0;
        for (int j = 0; j < length; ++j)
          if(max<m[i*length+j])
            max = m[i*length+j];
        for (int j = 0; j < length; ++j)
        {
            m[i*length+j] = max;
        }

    }
    MPI_Barrier(comm);
    MPI_Comm_split(comm,row,rank,&c1);
    MPI_Allreduce(m.data(),p_prime.data(),(length*length),MPI_INT,MPI_MAX,c1);
    MPI_Barrier(comm);
    sleep(rank);

    std::cout<<"PPrime : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << p_prime[i * length + j] << " ";
        std::cout << std::endl;
    }


    // Tree Hanging step

    sleep(rank);

    MPI_Barrier(comm);
    // for(int i=0;i<length;i++)
    // {
    //   for(int j=0;j<length;j++)
    //   {
    //     if(p_prime[j*length+i]==i+length*row)
    //     {
    //       m[j*length+i] = p[j*length+i];
    //     }
    //     else
    //     {
    //       m[j*length+i] = 0;
    //     }
    //   }
    // }

    for(int i = 0; i < length ; i++)
    {
        for(int j = 0; j < length; j++)
        {
            if(p[j+(length * i)] == i + ((rank / q) * length))
                m[j+(length * i)] = p_prime[j+(length * i)];
            else 
                m[j+(length * i)] = 0;
            if((i + ((rank / q) * length)) == 6)
                std::cout<<"P: "<<p[j+(length * i)]<<"::"<<"PPrime: "<<p_prime[j+(length * i)]<<std::endl;
        }
    }

    MPI_Barrier(comm);
    sleep(rank);

    std::cout<<"MTreeHanging : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << m[i * length + j] << " ";
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
    sleep(rank);

    for (int i = 0; i < length; ++i) {
        max=0;
        for (int j = 0; j < length; ++j)
          if(max<m[i*length+j])
            max = m[i*length+j];
        for (int j = 0; j < length; ++j)
        {
            m[i*length+j] = max;
        }
    }

    MPI_Barrier(comm);
    MPI_Comm_split(comm,row,rank,&c1);
    MPI_Allreduce(m.data(),c.data(),(length*length),MPI_INT,MPI_MAX,c1);
    MPI_Barrier(comm);
    sleep(rank);

    std::cout<<"CTreeHanging : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << c[i * length + j] << " ";
        std::cout << std::endl;
    }

    //std::cout<<"Hello1\n";
    MPI_Barrier(comm);
    sleep(rank);
    //std::cout<<"Hello2\n";
    for(int i=0 ; i< length; i++)
    {
        for(int j = 0; j< length; j++)
        {
        //std::cout<<"Hello3\n";
        if(p_prime[i*length + j] < c[i*length + j])
            pTreeHanging.push_back(c[i*length + j]);
        else 
            pTreeHanging.push_back(p_prime[i*length + j]);
        }
    }
    //std::cout<<"Hello4\n";
    MPI_Barrier(comm);
    sleep(rank);

    std::cout<<"PTreeHanging : \n";
    std::cout<<"Rank :"<<rank<<"\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j)
        std::cout << pTreeHanging[i * length + j] << " ";
        std::cout << std::endl;
    }

    return -1;
} // connected_components

#endif // A1_HPP


//ubuntu@ubuntu-xenial:/vagrant/A1/MPI_Testing$ mpic++ MPI_ScatterGather.cpp -std=c++11 -o MPI_ScatterGather
//ubuntu@ubuntu-xenial:/vagrant/A1/MPI_Testing$ mpirun -n 4 -mca btl ^openib ./MPI_ScatterGather 