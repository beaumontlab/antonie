// FOGSAA algorithm implementation with affine gap penalty
//               Developed  By Angana Chakraborty



// -------------------------------------------
// ############      DNA/Gene Sequences    ###############
/*-----------------------------------------
  Command line parameters (for alignments without affine gap penalty)
  ./fogsaa f1 f2 d 0 mt mst gp
  ---- where f1 --> name of the 1st file.
  f2 --> name of the 2nd file
  d  --> 1 for DNA sequences 
  0  --> as without affine gap penalty  
  mt --> score of a match
  mst--> score for mismatch
  gp --> gap penalty
  example ./fogsaa seq1.txt seq2.txt 1 0 1 -1 -2
  -----------------------------------------------------------
  -----------------------------------------------------------
  -----------------------------------------------------------
  Command line parameters ( for alignments with affine gap penalty)
  ./fogsaa f1 f2 d 1 mt mst go ge 
  ---- where f1 --> name of the 1st file.
  f2 --> name of the 2nd file
  d  --> 1 for DNA sequences 
  1  --> as with affine gap penalty 
  mt --> score of a match
  mst--> score for mismatch
  go --> gap open penalty
  ge --> gap extension penalty, Total penalty for a gap of length L= (go+L*ge)
  example ./fogsaa seq1.txt seq2.txt 1 1 1 -1 -10 -1
  --------------------------------------------------------------------
  --------------------------------------------------------------------

*/


#include<iostream>
#include<stdio.h>
#include<string.h> 
#include<stdlib.h>
using namespace std;
#include <sys/time.h>
#include <unistd.h>


int expanded=0,score=0,lowerBound,a[3]={0,0,0},b[3]={0,0,0},ch[3],optp1,optp2,approximate=0,threshold, affine=0, Mt,Mst,Gp,Go,Ge; // a=lower b=upper c=child index
// Mt == match score, Mst== Mismatch score, Gp== Gap penalty,, Go== Gap open penalty,, Ge== Gap extension penalty, 

// structure definition
int fscore[20][20],prmax,prmin;
char amino[25];
typedef struct celltype
{
  int present_score,lower,upper,type,filled; // type '1' is for (i+1,j+1),, type '2' is for (i,j+1), type '4' is for (i+1,j)  to identify them distinctly
}cell;

cell **c;


typedef struct qtype
{
  int p1,p2,type_upto_next,next_type,next_lower;
  struct qtype *next;
}queue;
int maxPointer=-1;
queue **q=NULL;

string out1, out2;
string pair1, pair2;
unsigned int curp1=0,curp2=0;

// function prototype
int ins_queue(int,int,int ,int,int ,int);
void del_queue(void);
void align(unsigned int, unsigned int);
void sort(void);
void calculate_score(int p,int *l,int *u,int p1,int p2);

int main(int argc, char * argv[])
{
  char ca,cb,pathend=0;
  int type_total=1,new_type=0,new_score=0,np1=0,np2=0,new_lower,new_upper,next_lower,next_upper;
  int p,present;

  //            1   2  3
  // DNA mode: in1 in2 1 affine-bit match-score mismatch-score gap-penalty
  // parameter checking

  affine = argv[4]==string("1");
  
  if(!affine)
    {  
      if(argc!=8)
	printf(" Invalid argument\n");
      else
	{
	  Mt=atoi(argv[5]);
	  Mst=atoi(argv[6]);
	  Gp=atoi(argv[7]);
	}
    }
  else if(affine) // gene seq with affine gap
    {
      if(argc!=9)
	printf(" Invalid argument\n");
      else
	{
	  Mt=atoi(argv[5]);
	  Mst=atoi(argv[6]);
	  Go=atoi(argv[7]);
	  Ge=atoi(argv[8]);
	}
    }

  pair1=argv[1];
  pair2=argv[2];
  unsigned int i;

  c = (cell **)malloc((pair1.size()+1) * sizeof(cell *));
  for (i = 0; i <= pair1.size(); i++) {
    c[i] =(cell *) calloc( (pair2.size()+1), sizeof(cell));
  }
	
  c[0][0].present_score=0;
  c[0][0].type=-1;
	
  calculate_score(c[0][0].present_score,&(c[0][0].lower),&(c[0][0].upper),0,0);
  lowerBound=c[0][0].lower;
	

  // estimate the threshold value with 30% similarity requirement
  int ml,sl,th;
  if(pair1.size() > pair2.size())
    {
      ml=pair1.size();
      sl=pair2.size();
    }
  else
    {
      ml=pair2.size();
      sl=pair1.size();
    }
  th=ml*30/100;
  if(affine==0)
    threshold=th*Mt+(sl-th)*Mst+Gp*(ml-sl);
  else
    threshold=th*Mt+(sl-th)*Mst+Ge*(ml-sl)+Go;
      
  unsigned int v=c[0][0].upper-c[0][0].lower+1;
  cout<<"Creating queue of size "<<v<<endl;
  q=new queue*[10*v];  // creating array of pointers to queuetype
  for(i=0;i<10*v;i++)
    q[i]=NULL;

  if(pair1.size()!=0 && pair2.size()!=0)
    {
      do
	{
	  pathend=1;
	  while((curp1<=pair1.size()-1) || (curp2 <=pair2.size()-1))
	    {
	      present=c[curp1][curp2].present_score;
	      // expand the child node
			
	      if((type_total==1)||(type_total==2)||(type_total==4))
		{ 
		  //this is the first child
				
		  if(curp1<=pair1.size()-1 && curp2<=pair2.size()-1)
		    {
					
		      ca=pair1[curp1];
		      cb=pair2[curp2];
		      if(ca==cb)
			p=Mt;
		      else
			p=Mst;

			                 
					
		      calculate_score(present+p,&a[0],&b[0],curp1+1,curp2+1);
		      if(affine==0)
			{
					
			  calculate_score(present+Gp,&a[1],&b[1],curp1,curp2+1);
			  calculate_score(present+Gp,&a[2],&b[2],curp1+1,curp2);
			}
		      else // affine gap penalty
			{
			  if((c[curp1][curp2].type==1)||(c[curp1][curp2].type==-1))
			    {
			      calculate_score(present+Go+Ge,&a[1],&b[1],curp1,curp2+1);
			      calculate_score(present+Go+Ge,&a[2],&b[2],curp1+1,curp2);
			    }
			  else if(c[curp1][curp2].type==2)// as the gap is already opened in the first chain
			    {
			      calculate_score(present+Ge,&a[1],&b[1],curp1,curp2+1);
			      calculate_score(present+Go+Ge,&a[2],&b[2],curp1+1,curp2);
			    }
			  else // as the gap is already opened in the 2nd chain
			    {
			      calculate_score(present+Go+Ge,&a[1],&b[1],curp1,curp2+1);
			      calculate_score(present+Ge,&a[2],&b[2],curp1+1,curp2);
			    }
						
			}
					
					
		      sort(); // sort the children according to upper and when there is tie... check the lower
					
					
					
		      if(ch[0]==1)
			{
						
			  // child would be of type 1
			  new_type=1;
			  np1=curp1+1;
			  np2=curp2+1;
			  new_score=present+p;
						
						

			}
		      else if(ch[0]==2)
			{	
						
			  // child would be of type 2
			  new_type=2;
			  np1=curp1;
			  np2=curp2+1;
			  if(affine==0)
			    new_score=present+Gp;
			  else // in the affine gap penalty scheme
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 ||(c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
				new_score=present+Go+Ge;
			      else if(c[curp1][curp2].type==2)// the gap was already opened in first chain
				new_score=present+Ge;
							
							
			    }
						
						


						
			}
		      else
			{
			  // child would be of type 4
						
			  new_type=4;
			  np1=curp1+1;
			  np2=curp2;
			  if(affine==0)
			    new_score=present+Gp;
			  else // in the affine gap penalty scheme
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now in 1st chain
						
				new_score=present+Go+Ge;
			      else // the gap was already opened in 2nd chain
				new_score=present+Ge;
			    }
			}

		      //printf("before queue\n");
		      maxPointer=ins_queue(curp1,curp2,new_type+ch[1],ch[1],a[1],b[1]);
		      //printf("after queue\n");
		    }
		  else if(curp1<=pair1.size()-1)
		    {
		      // only type '4' child is possible
		      new_type=4;
		      np1=curp1+1;
		      np2=curp2;
		      if(affine==0)
			new_score=present+Gp;
		      else // in the affine gap penalty scheme
			{
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
			    new_score=present+Go+Ge;
			  else // the gap was already opened
			    new_score=present+Ge;
							
							
			}
						
				
		    }
		  else
		    {
		      // only type '2' child is possible
		      new_type=2;
		      np1=curp1;
		      np2=curp2+1;
		      if(affine==0)
			new_score=present+Gp;
		      else // in the affine gap penalty scheme
			{
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
			    new_score=present+Go+Ge;
			  else // the gap was already opened
			    new_score=present+Ge;
							
							
			}
						
				

		    }
			
			
		}  // end of first child
	      else if((type_total==3)||(type_total==5)||(type_total==6))
		{
		  // this is the 2nd child
				
		  if(new_type==1)
		    {
		      // 2nd child is of type type 1
		      np1=curp1+1;
		      np2=curp2+1;
		      ca=pair1[curp1];
		      cb=pair2[curp2];
		      if(ca==cb)
			p=Mt;
		      else
			p=Mst;
		      new_score=present+p;
		      if(7-type_total==2)
			{
			  if(affine==0)
			    calculate_score(present+Gp,&next_lower,&next_upper,curp1,curp2+1);
			  else
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
				calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1,curp2+1);
			      else
				calculate_score(present+Ge,&next_lower,&next_upper,curp1,curp2+1);
			    }
			}
		      else if(7-type_total==4)
			{
			  if(affine==0)
			    calculate_score(present+Gp,&next_lower,&next_upper,curp1+1,curp2);
			  else
						
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
				calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1+1,curp2);
			      else
				calculate_score(present+Ge,&next_lower,&next_upper,curp1+1,curp2);
							
			    }
			}
		      //printf("before queue2\n");
		      maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
		      //printf("before queue2\n");

		    }
		  else if(new_type==2)
		    {
		      //2nd child is of type type 2
		      np1=curp1;
		      np2=curp2+1;
		      if(affine==0)
			new_score=present+Gp;
		      else
			{	
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
			    new_score=present+Go+Ge;
			  else 
			    new_score=present+Ge;  // because gap has already been opened
							
			}
		      if(7-type_total==1)
			{
			  ca=pair1[curp1];
			  cb=pair2[curp2];

			  if(ca==cb)
			    p=Mt;
			  else
			    p=Mst;
			  calculate_score(present+p,&next_lower,&next_upper,curp1+1,curp2+1);
			}
		      else if(7-type_total==4)
			{
			  if(affine==0)
			    calculate_score(present+Gp,&next_lower,&next_upper,curp1+1,curp2);
			  else
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
				calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1+1,curp2);
			      else
				calculate_score(present+Ge,&next_lower,&next_upper,curp1+1,curp2);
			    }
			}
		      //printf("before queue3\n");
		      maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
		      //printf("before queue3\n");
					
		    }
		  else
		    {
		      // 2nd child is of type type 4
		      np1=curp1+1;
		      np2=curp2;
		      if(affine==0)
			new_score=present+Gp;
		      else
			{	
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
			    new_score=present+Go+Ge;
			  else 
			    new_score=present+Ge;  // because gap has already been opened
							
			}
		      if(7-type_total==1)
			{
			  ca=pair1[curp1];
			  cb=pair2[curp2];
			  if(ca==cb)
			    p=Mt;
			  else
			    p=Mst;
			  calculate_score(present+p,&next_lower,&next_upper,curp1+1,curp2+1);
			}
		      else if(7-type_total==2)
			{
			  if(affine==0)
			    calculate_score(present+Gp,&next_lower,&next_upper,curp1,curp2+1);
			  else
			    {
			      if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
				calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1,curp2+1);
			      else
				calculate_score(present+Ge,&next_lower,&next_upper,curp1,curp2+1);
			    }
			}
		      //printf("before queue4\n");
		      maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
		      //printf("before queue4\n");

		    }



		}// end of 2nd child
	      else if(type_total==7)
		{
		  // this is the third child
				

		  if(new_type==1)
		    {
		      np1=curp1+1;
		      np2=curp2+1;
		      ca=pair1[curp1];
		      cb=pair2[curp2];
		      if(ca==cb)
			p=Mt;
		      else
			p=Mst;
		      new_score=present+p;
		    }
		  else if(new_type==2)
		    {
		      np1=curp1;
		      np2=curp2+1;
		      if(affine==0)
			new_score=present+Gp;
		      else
			{	
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
			    new_score=present+Go+Ge;
			  else 
			    new_score=present+Ge;  // because gap has already been opened
							
			}
		    }
		  else
		    {
		      np1=curp1+1;
		      np2=curp2;
		      if(affine==0)
			new_score=present+Gp;
		      else
			{	
			  if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
			    new_score=present+Go+Ge;
			  else 
			    new_score=present+Ge;  // because gap has already been opened
							
			}
		    }

		  // here no node is to be inserted in the queue
		}
			

			
	      // write the new child node
	      cout<<np1<<" "<<np2<<endl;
	      cout<<c[np1][np2].type<<endl;
	      cout<<c[np1][np2].present_score<<endl;
	      cout<<c[np1][np2].filled<<endl;
	      cout<<new_score<<endl;
	      
	      if(c[np1][np2].type<=4 && c[np1][np2].present_score>=new_score && c[np1][np2].filled==1) // skip the path if already expanded (same)node is better
		{
		  // printf("pruned already better at(%d,%d) old=%d,new=%d\n",np1,np2,c[np1][np2].present_score,new_score);
		  pathend=0;
					
		  break;
		}
		
	      else  // insert in the cell
		{
		  calculate_score(new_score,&new_lower,&new_upper,np1,np2);
		  c[np1][np2].present_score=new_score;
		  c[np1][np2].lower=new_lower;
		  c[np1][np2].upper=new_upper;
		  c[np1][np2].type=new_type;
		  c[np1][np2].filled=1;
				
		}

		
	      // print the node
				
		
	      /*char c1,c2;
			
		ca=pair1[np1-1];
		cb=pair2[np2-1];
		if(np1==0)
		c1=' ';
		else if(curp1==np1)
		c1='-';
		else
		c1=ca;

		if(np2==0)
		c2=' ';
		else if(curp2==np2)
		c2='-';
		else
		c2=cb;
		printf("node(%c,%c,%d,%d,%d,%d,%d)\n",c1,c2,np1,np2,c[np1][np2].present_score,c[np1][np2].lower,c[np1][np2].upper);*/
			
			
	      // now point the child as the current node
	      curp1=np1;
	      curp2=np2;
	      type_total=1;

	      if(c[np1][np2].upper<lowerBound)  
		{
		  pathend=0;
		  break;
		}
			
			
	      expanded++; // counts no of nodes expanded


			
      

	    }// end of while

		

	  if((c[curp1][curp2].present_score>lowerBound) && pathend==1) // if the current path is not totally expanded , don't reset the lowerbound
	    {
	      lowerBound=c[curp1][curp2].present_score;
	      optp1=curp1;
	      optp2=curp2;
	      //printf("Path end : new lowerBound=%d\n",lowerBound);
	    }
	
	  /*if((expanded>m*n/3) && (lowerBound<threshold))
	    {
	    printf("The given sequences are not globally similar(below 30%% of similarity)-- try local allignment --\n");
	    printf(" The score below is near optimal\n");
	    printf("The threshold== %d\n",threshold);
	    approximate=1;
	    printf("now the future is= %d\n",new_upper);
	    break;
	    }*/


	
	
	  if(maxPointer!=-1 )
	    {
	      curp1=q[maxPointer]->p1;
	      curp2=q[maxPointer]->p2;
				
	      type_total=q[maxPointer]->type_upto_next;
	      new_lower=q[maxPointer]->next_lower;
	      new_upper=maxPointer+c[0][0].lower;
	      new_type=q[maxPointer]->next_type;
	      del_queue();
				
				
	    }
	  //printf(" The future is=%d\n",new_upper);
	  int cpl;
	  if(curp1>curp2)
			
	    cpl=curp1;
				
	  else
	    cpl=curp2;
				
	  if(((cpl>70*ml/100)&&(lowerBound<threshold))||(new_upper<threshold))
	    {
	      printf("The given sequences are not globally similar(below 30%% of similarity)-- try local allignment --\n");
	      printf(" The score below is near optimal\n");
	      printf("The threshold== %d\n",threshold);
	      approximate=1;
	      printf("now the future is= %d\n",new_upper);
	      break;
	    }
		
	}while(lowerBound<new_upper);



      align(optp1,optp2);
      cout<<out1<<endl;
      cout<<out2<<endl;
      
      printf("total nodes expanded==%d\n\n",expanded);
			
      printf("score= %d\n",c[optp1][optp2].present_score);
    }// end of if
}// end of main



void align(unsigned int p1, unsigned int p2)
{
  char c1,c2;
  unsigned int pp1,pp2;

  if(c[p1][p2].type==-1)
    {
      return;
    }
  else
    {
      if(c[p1][p2].type==1)
	{
	  pp1=p1-1;
	  pp2=p2-1;
	}
			
      else if(c[p1][p2].type==2)
	{
	  pp1=p1;
	  pp2=p2-1;
	}
      else
	{
	  pp1=p1-1;
	  pp2=p2;
	}

      align(pp1,pp2);
      if(p1==pair1.size()+1)
	{
	  out1.append(1,'-');
	  c1='-';
	}
      else
	{
	  if(pp1==p1)
	    {
	      out1.append(1,'-');	      
	      c1='-';
	    }
	  else
	    {
	      char ca=pair1[p1-1];
	      out1.append(1, ca);
	      c1=ca;
	    }
	}

      if(p2==pair2.size()+1)
	{
	  out2.append(1,'-');
	  c2='-';
	}
      else
	{
	  if(pp2==p2)
	    {
	      out2.append(1,'-');
	      c2='-';
	    }
	  else
	    {
	      char cb=pair2[p2-1];
	      out2.append(1,cb);
	      c2=cb;
	    }
	}
		
      
    } // end of else

}


void del_queue()
{
  queue * temp;
  if(maxPointer==-1)
    return;
  else
    {
      temp=q[maxPointer];
      q[maxPointer]=temp->next;
      free(temp);
      while(q[maxPointer]==NULL)  // this row is empty , point next row
	maxPointer--;
    }
}


int ins_queue(int p1,int p2,int type_total,int next_type,int next_lower,int next_upper)
{
  int inserted=0;
	
  queue *p,*prev,*newnode;
  if((next_upper-c[0][0].lower)>=0)
    {
      newnode=(queue*)calloc(1, sizeof(queue));
      newnode->p1=p1;
      newnode->next_type=next_type;
      newnode->next_lower=next_lower;
      newnode->type_upto_next=type_total;
      newnode->p2=p2;
      printf("queue insertion of =(%d,%d) type %d at %d\n",p1,p2,next_type,next_upper);
	
      if(maxPointer==-1)
	{
	  q[next_upper-c[0][0].lower]=newnode;
	  newnode->next=NULL;
	  maxPointer=next_upper-c[0][0].lower;
	}
      else
	{
	  if(q[next_upper-c[0][0].lower]==NULL)  // first memeber in this row
	    {
	      q[next_upper-c[0][0].lower]=newnode;
	      newnode->next=NULL;
	      if((next_upper-c[0][0].lower)>maxPointer)
		maxPointer=next_upper-c[0][0].lower;
	    }
	  else
	    {
	      // search the appropriate position in the row comparing the lower value
	      cout<<"Reading q["<<next_upper-c[0][0].lower<<"]"<<endl;
	      p=q[next_upper-c[0][0].lower];
			
	      prev=NULL;
	      while(p!=NULL)
		{
		  if(p->next_lower<=next_lower)
		    {
		      // insert here(before p)
		      if(prev==NULL) // insert as the first row in the list
			{
			  q[next_upper-c[0][0].lower]=newnode;
			  newnode->next=p;
			}
		      else
			{
			  // insert in the middle of the row
			  prev->next=newnode;
			  newnode->next=p;
			}
		      inserted=1;
		      break;


		    }
		  else
		    {
		      prev=p;
		      p=p->next;
		    }
		}

	      if(inserted==0)
		{
		  // insert at the end
		  prev->next=newnode;
		  newnode->next=NULL;
		}
	    }
	

	}
    }
  //else printf("skipped************************\n");
  return(maxPointer);
	
} // end of the function

void calculate_score(int p,int* l,int * u,int p1,int p2)
{
  if(affine==1)
    {
      Gp=Go+Ge;
    }
  if((pair1.size()-p1)<=(pair2.size()-p2))
    {
      *l=p+((pair1.size()-p1)*Mst+Gp*((pair2.size()-p2)-(pair1.size()-p1)));
      *u=p+((pair1.size()-p1)*Mt+Gp*((pair2.size()-p2)-(pair1.size()-p1)));
    }
  else
    {
      *l=p+Mst*(pair2.size()-p2)+Gp*((pair1.size()-p1)-(pair2.size()-p2));
      *u=p+(pair2.size()-p2)*Mt+Gp*((pair1.size()-p1)-(pair2.size()-p2));
    }
}

void sort(void)  // descending order
{
  int i,j,t;
  ch[0]=1;
  ch[1]=2;
  ch[2]=4;
  for(i=0;i<2;i++)
    {
      for(j=0;j<2-i;j++)
	{
	  if((a[j]<a[j+1])||((a[j]==a[j+1])&&(b[j]<b[j+1])))
	    {
	      //swap
	      t=a[j];
	      a[j]=a[j+1];
	      a[j+1]=t;

	      t=ch[j];
	      ch[j]=ch[j+1];
	      ch[j+1]=t;

	      t=b[j];
	      b[j]=b[j+1];
	      b[j+1]=t;

	    }
			
	}
    }
}






