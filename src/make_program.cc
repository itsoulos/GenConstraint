# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>

const char *short_options="hp:o:";
char problemfile[1024];
char outputfile[1024];

void	print_usage()
{
       printf("\t-h    	Print this help screen.\n"
		"\t-p	Specify problem file.\n"
		"\t-o	Specify output file.\n");
}

void	parse_cmd_line(int argc,char **argv)
{
	if(argc==1)
	{
		print_usage();
		exit(EXIT_FAILURE);
	}
	int next_option;
	do
	{
		next_option=getopt(argc,argv,short_options);
		switch(next_option)
		{
			case 'p':
				strcpy(problemfile,optarg);
				break;
			case 'o':
				strcpy(outputfile,optarg);
				break;
			case -1:
				break;
			case '?':
				print_usage();
				exit(EXIT_FAILURE);
				break;
			default:
				print_usage();
				exit(EXIT_FAILURE);
				break;
		}
	}while(next_option!=-1);
}

int main(int argc,char **argv)
{
	strcpy(outputfile,"constraint");
	strcpy(problemfile,"");
	parse_cmd_line(argc,argv);
	if(strlen(problemfile)==0) return 1;
	int pos_of_dot=0;
	
	for(int i=0;i<strlen(problemfile);i++)
	{
		if(problemfile[i]=='.') 
		{
			pos_of_dot=i;
			break;
		}
	}
	char extension[100];
	char compiler[100];
	strcpy(extension,&problemfile[pos_of_dot+1]);
	FILE *fp;
	if(!strcmp(extension,"cc") || 
			!strcmp(extension,"c++")
			||!strcmp(extension,"CC")
			||!strcmp(extension,"cpp"))
	{
		strcpy(compiler,CXX);
	}
	else
	if(!strcmp(extension,"f") ||
		!strcmp(extension,"F")
		||!strcmp(extension,"for"))
	{
		strcpy(compiler,F77);
	}
	else
	if(!strcmp(extension,"c"))
	{
		strcpy(compiler,CC);
	}
	char compile_command1[1024];
	char compile_command2[1024];
	char rootdir[1024];
	strcpy(rootdir,ROOTDIR);

	char	files[7][1024]={"tolmin.o","base.o","doublepop.o",
				"main.o","problem.o","getoptions.o","cobyla.o"};
	char	cmdline[2048];
	strcpy(cmdline,"");
	for(int i=0;i<7;i++)
	{
		sprintf(cmdline,"%s %s/src/%s ",cmdline,rootdir,files[i]);
	}

	unlink("myproblem.o");
	if(!strcmp(compiler,F77))
		sprintf(compile_command1,"%s %s %s -c %s -o myproblem.o",compiler,
				F77FLAGS,INCL,problemfile);
	else
		sprintf(compile_command1,"%s %s  -c %s -o myproblem.o",compiler,
				INCL,problemfile);


	int code1=system(compile_command1);
	if(code1==-1) 
	{
		return printf("THE COMPILATION OF %s HAS BEEN FAILED\n",problemfile);
	}

	sprintf(compile_command2,"%s %s myproblem.o %s -o %s",CXX,cmdline,LINK,outputfile);

	unlink(outputfile);

	int code2=system(compile_command2);
	if(code2==-1) 
	{
		return printf("THE CREATION OF %s HAS BEEN FAILED\n",outputfile);
	}

	FILE *test=fopen(outputfile,"rb");
	if(test==NULL) 
		return printf("EXECUTABLE FILE WAS NOT CREATED\n");
	fclose(test);

	printf("RUN ./%s IN ORDER TO RUN THE PROBLEM\n",outputfile);
	return 0;
}
