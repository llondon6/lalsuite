/*
   Parse Command line:
  
   Typical Command line:

   laldoc inputFile.c errorFile /home/alice/errorDir/ /home/alice/inputDir/
     -0-     -1-          -2-           -3- (optnl)         -4- (optnl)
  
   Obviously 0,1 and 2 are necessary: The executable, the source input
   and the error file name.  

   Number 3 is the directory where the the error file will appear.
   If this is not specified, the errorFile will be opened the current
   directory.

   Number 4 is the directory where the program will look for the 
   source input. If this is not specified, the program will look
   in the current directory.


*/

#include "LALDoc.h"

int ParseCommandLine( int  argcount , char **argvector , char *source , 
                      FILE **ptrptrInputFile )
{
        int i,dirLen;
        char ErrFileName[MAXSTR],ErrDir[MAXSTR],FullErrStr[MAXSTR];
        char SourceDir[MAXSTR],FullSourceStr[MAXSTR];
        for (i=0;i<MAXSTR;i++){
                ErrFileName[i]      = '\0';
                ErrDir[i]           = '\0';
                FullErrStr[i]       = '\0';
                SourceDir[i]        = '\0';
                FullSourceStr[i]    = '\0';
        }

        if( strcmp(argvector[1], "-h") == 0 ){
                fprintf(stderr ,"Help: Must specify an input and exception file\n");
                exit(BADEXIT);
        }

        if(argcount < 3 ){
                fprintf(stderr , "Need both source input and error.\n");
                exit(BADEXIT);
        }

        strcpy( source        , argvector[1] ) ;
        strcpy( ErrFileName   , argvector[2] );


        /* 
         * Opening the error reporting file 
        */
        /* If there is a directory specified for the error file ... */
        dirLen = 0;
        if (argcount >= 4 ){
                strcpy( ErrDir , argvector[3] ) ;
                dirLen = strlen(ErrDir);
                i=0;
                /* ...  copy it into a string  and ... */
                while(*(ErrDir+i) != '\0' ){
                        *(FullErrStr+i) = *(ErrDir+i);
                        i++;
                }
        }
        /* ... copy the error file name onto the end of the directory and  ... */
        i=0;
        while( *(ErrFileName+i)!= '\0' ){
               *(FullErrStr+dirLen+i)=*(ErrFileName+i);
               i++;
        }
        /* ... open the error reporting file. */
        ptrLALErrorFile = OpenAFile( FullErrStr , "a" , 1 );
        fprintf(ptrLALErrorFile,"%%Laldoc: Reporting Errors in: %s\n",source);        



        /* 
         * Opening the input file 
        */
        /* If there is a directory specified for the input file ... */
        dirLen = 0;
        if (argcount >= 5 ){
                strcpy( SourceDir , argvector[4] ) ;
                dirLen = strlen(SourceDir);
                i=0;
                /* ...  copy it into a string  and ... */
                while(*(SourceDir+i) != '\0' ){
                        *(FullSourceStr+i) = *(SourceDir+i);
                        i++;
                }
        }
        /* ... copy the input file name onto the end of the directory and  ... */
        i=0;
        while( *(source+i)!= '\0' ){
               *(FullSourceStr+dirLen+i)=*(source+i);
               i++;
        }
        /* ... open the input reporting file. */
        *ptrptrInputFile = OpenAFile( FullSourceStr , "r" , 0 );

        if( (*ptrptrInputFile)==NULL ) {
               fprintf(stderr,"Couldn't find Source file: %s\n",FullSourceStr);
               exit(BADEXIT);
                
        }

        return 1 ;
}
