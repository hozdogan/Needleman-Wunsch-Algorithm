/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package paralel;

 import java.lang.*;
import java.util.*;
import java.io.*;

public class Paralel {

   
    public static void main(String[] args) {
        float [][] A=new float[202][202];
        String [] sekans5k = new String[5000];
        String [] sekans10k = new String[10000];
        String [] sekans15k = new String[15000];
        float [][] skors = new float[5000][5000];
        float [][] tablo = new float[20][3];
        //char [] nucleotits={'A','T','G','C'};
       
        //String dynmc="";
        
       // String sequence1 = "ATCGTACGTACGTTGCATGCAGGTACTTAACCTGACTGAACTGACTGGTACTGCCATGAACCCAAAGGTCCAGTATCGATCGATCGGGTCGGGCTGGGACCAGTGCTGAAAACCTGGCTGAGCCATGGCTGGAGAATTTCGGATGCGGTAGGAGTCGATAGCTGAGTCGTAGCTGAGTGAGTCGATGCTAGGTAGTTGAT";
        //String sequence2 = "GGTAGCTGATGGGCTAGTTTATATGCTTGGATTGCTATGGGATGCTTTGATTGCGTAGTTTCGATTCGGTAGCTGGGGTAAAATTTCGGCTGAGCTGGATCGGCTCGATTGCTAGGTCGATGGCTAGTGCTAGGTCGTAGTCGTTTAGCTGGAATCGATGGATCGGGGGATGCGTAGTTGGGCTAGTTCGATGCTAGTCG";
        
       /* for(int i=0;i<sekans15k.length;i++)//5000 10000 15000 değiştir
        {
            for(int j=0;j<200;j++)
            {
                Random r = new Random();
                int a=r.nextInt(4);
                dynmc+=nucleotits[a];
               
            }
            sekans15k[i]=dynmc;
            dynmc="";
        }*/
        /*for(int k=0;k<5000;k++)
        {
            System.out.println(sekans5k[k]);
            
        }*/
        
        int index5k=0,index10k=0,index15k=0;
        try
        {
             File out = new File("C:\\Users\\asus\\Desktop\\5K_Sequence.fasta");
             Scanner sc = new Scanner(out);
                while (sc.hasNextLine()) {
                String line = sc.nextLine();
                
                if(line.startsWith(">"))
                {
                    line =sc.nextLine();
                    sekans5k[index5k]=line;
                    index5k++;
                   
                }
                
                
              }
        }catch(Exception e){e.printStackTrace();}
        
        
        try
        {
             File out = new File("C:\\Users\\asus\\Desktop\\10K_Sequence.fasta");
             Scanner sc = new Scanner(out);
                while (sc.hasNextLine()) {
                String line = sc.nextLine();
                
                if(line.startsWith(">"))
                {
                    line =sc.nextLine();
                    sekans10k[index10k]=line;
                    index10k++;
                   
                }
                
                
              }
        }catch(Exception e){e.printStackTrace();}
        
        
        
        try
        {
             File out = new File("C:\\Users\\asus\\Desktop\\15K_Sequence.fasta");
             Scanner sc = new Scanner(out);
                while (sc.hasNextLine()) {
                String line = sc.nextLine();
                
                if(line.startsWith(">"))
                {
                    line =sc.nextLine();
                    sekans15k[index15k]=line;
                    index15k++;
                   
                }
                
                
              }
        }catch(Exception e){e.printStackTrace();}
        
        
        /*for(int k=0;k<5000;k++)
        {
            System.out.println(sekans5k[k]);
            
        }
        */
        //String sequence1="ATGCAT";
        //String sequence2="GTGCAT";
        float gap=-1.832482334f,match=3.621354295f,mismatch=-2.451795405f,sol,üst,capraz,max,skor,elem;
       
      
        for(int i=0;i<A.length;i++)
        {
            for(int j=0;j<A[i].length;j++)
            {
                A[i][j]=0;
            }
        }
        int l=0;
        A[0][1]='-';//-
        A[1][0]='-';//-
       
        //System.out.println(sequence1.charAt(l));
        
       /* for(int i=2;i<A.length;i++)//sekansları yerleştirme kısmı
        {
            A[i][0]=sequence1.charAt(l);
            A[0][i]=sequence2.charAt(l);
            l++;
        }       */
        
       for(int r=0;r<200;r++)//sekans5k.length
       {
           for(int c=r+1;c<200;c++)
           {
               
                for(int i=2;i<A.length;i++)//sekansları yerleştirme kısmı
                {
                    A[i][0]=sekans5k[r].charAt(l);
                    A[0][i]=sekans5k[c].charAt(l);
                l++;
                }     
           l=0;
           for(int row=1;row<A.length;row++)
            {
            for(int col=1;col<A[row].length;col++)
            {
                elem=A[row][col];
                
                if(A[0][col]=='-'||A[row][0]=='-')
                {
                    if(A[row][0]=='-')
                    {
                        elem=A[row][col-1]+gap;
                      
                    }
                    else if(A[0][col]=='-')
                    {
                        elem=A[row-1][col]+gap;
                        
                    }
                }
                if(A[0][col]=='-'&&A[row][0]=='-')
                {
                    elem=0;
                }
                A[row][col]=elem;
            }
        }
        for(int row=2;row<A.length;row++)
        {
            for(int col=2;col<A[row].length;col++)
            {
                elem=A[row][col];
                if(A[0][col]!=A[row][0])
                {
                    sol=A[row][col-1]+gap;
                    üst=A[row-1][col]+gap;
                    capraz=A[row-1][col-1]+mismatch;
                    max=sol;
                    if(üst>max)
                    {
                        max=üst;
                    }
                    if(capraz>max)
                    {
                        max=capraz;
                    }
                    elem=max;
                    
                }
                else if(A[0][col]==A[row][0])
                {
                    sol=A[row][col-1]+gap;
                    üst=A[row-1][col]+gap;
                    capraz=A[row-1][col-1]+match;
                    max=sol;
                    if(üst>max)
                    {
                        max=üst;
                    }
                    if(capraz>max)
                    {
                        max=capraz;
                    }
                    elem=max;
                    
                }
                A[row][col]=elem;
            }
        }
               
               skor=A[A.length-1][A.length-1];
             
               skors[r][c]=skor;
              System.out.println("Hizalama "+"x = "+r+" y = "+c+" Skor = "+skors[r][c]);
           }
         //System.out.println("Hizalama "+"x = "+r);
       }
       
       for(int x=0;x<20;x++)
       {
           Random rx = new Random();
           Random cx = new Random();
           int r=rx.nextInt(200);
           int c = cx.nextInt(200);
           while(r==c)
           {
                r=rx.nextInt(200);
                c = cx.nextInt(200);
           }
           tablo[x][0]=r;
           tablo[x][1]=c;
           if(skors[r][c]==0)
           {
                 tablo[x][2]=skors[c][r];
           }
           else if(skors[r][c]!=0)
           {
               tablo[x][2]=skors[r][c];
           }
          
           
       }
       for(int i=0;i<20;i++)
       {
          
           System.out.print("x = "+tablo[i][0]+" y = "+tablo[i][1]+" Skor = "+tablo[i][2]);
           
           System.out.print("\n");
           System.out.print("\n");
       }
        
      
}
}
