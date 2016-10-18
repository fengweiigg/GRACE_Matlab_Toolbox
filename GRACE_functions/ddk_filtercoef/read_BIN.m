%function which reads in a binary file containing symmetric/full or block
%diagonal matrices and associated vectors and parameters
%Roelof Rietbroek, 7-1-2008
%updated: 29-07-2008
%
%usage: dat=read_BIN(file)
%returns a structure array 'dat' with the file content
% the matrix remains in packed form (dat.pack1 field)
%or:    dat=read_BIN(file,'F')
% also expands the matrix to its full form (dat.mat1 field)
% Warning: this option may cause excessive RAM memory use with large matrices
%function now also works with Octave


function dat=read_BIN(file,varargin)

unpack=false;


for i=1:size(varargin,2)
  switch varargin{i}
   case {'F'}
    unpack=true; % unpack matrix in full size
  end
end



%open file for read acces in little endian format
[fid,message]=fopen(file,'r','ieee-le');
%check for errors
 if (fid == -1) 
   message
   dat=[];
   return; 
 end


%read BINARY version and type from file

dat.version(1,1:8)=fread(fid,8,'uint8=>char')';

dat.type(1,1:8)=fread(fid,8,'uint8=>char')';
dat.descr(1,1:80)=fread(fid,80,'uint8=>char')';

%read indices
%integers:inint,indbls,inval1,inval2,ipval1,ipval2
metaint=fread(fid,6,'integer*4');


%put index data in structure array
n_ints=metaint(1);
n_dbls=metaint(2);
dat.nval1=metaint(3);
dat.nval2=metaint(4);
dat.pval1=metaint(5);
dat.pval2=metaint(6);


%Type dependent index data
switch dat.type
 case {'BDSYMV0_','BDFULLV0'}
  %read additional nblocks parameter
  dat.nblocks=fread(fid,1,'integer*4');
end 

%get meta data 
%integers
if(n_ints > 0)
  list=fread(fid,n_ints*24,'uint8=>char');
  dat.ints_d=reshape(list,24,n_ints)';
  dat.ints=fread(fid,n_ints,'integer*4');
end

%doubles
if(n_dbls > 0)
  list=fread(fid,n_dbls*24,'uint8=>char');
  dat.dbls_d=reshape(list,24,n_dbls)';
  dat.dbls=fread(fid,n_dbls,'real*8');
end

%side description meta data
list=fread(fid,dat.nval1*24,'uint8=>char');
%reshape characters and put in dat struct array
dat.side1_d=reshape(list,24,dat.nval1)';



%type specific meta data

switch dat.type
 case {'BDSYMV0_','BDFULLV0'}
  %read additional nblocks parameter
  dat.blockind=fread(fid,dat.nblocks,'integer*4');
end 



%data (type dependent)

switch dat.type
  
 case 'SYMV0___'

  dat.pack1=fread(fid,dat.pval1,'real*8');

  if( unpack)
    row=[1:dat.nval1];
    col=row;
    %the array ind is upper triangular
    ind=triu(row'*ones(1,dat.nval1)) +triu(ones(dat.nval1,1)*(col.*(col-1))/2);
    %copy data from packed vector to full array
    
    
    dat.mat1=dat.pack1(ind+ind'-diag(diag(ind)));
    clear ind
    dat=rmfield(dat,'pack1');
  end
 case 'SYMV1___'
  dat.vec1=fread(fid,dat.nval1,'real*8');

  dat.pack1=fread(fid,dat.pval1,'real*8');

  if( unpack)
    row=[1:dat.nval1];
    col=row;
    %the array ind is upper triangular
    ind=triu(row'*ones(1,dat.nval1)) +triu(ones(dat.nval1,1)*(col.*(col-1))/2);
    %copy data from packed vector to full array
    
    dat.mat1=dat.pack1(ind+ind'-diag(diag(ind)));
    
    clear ind
    dat=rmfield(dat,'pack1');
  end
  
 case 'SYMV2___'
    dat.vec1=fread(fid,dat.nval1,'real*8');
    dat.vec2=fread(fid,dat.nval1,'real*8');
    dat.pack1=fread(fid,dat.pval1,'real*8');


  if(unpack)
    row=[1:dat.nval1];
    col=row;
    %the array ind is upper triangular
    ind=triu(row'*ones(1,dat.nval1)) +triu(ones(dat.nval1,1)*(col.*(col-1))/2);
    %copy data from packed vector to full array
    
    dat.mat1=dat.pack1(ind+ind'-diag(diag(ind)));
    
    clear ind
      dat=rmfield(dat,'pack1');
  end
    
 case 'BDSYMV0_'
    dat.pack1=fread(fid,dat.pval1,'real*8');
    
    if(unpack)
      dat.mat1=[];
      skip=0;
      skipentries=0;
      for i=1:dat.nblocks
        
        sz=dat.blockind(i)-skip;
        
        row=[1:sz];
        col=row;
        ind=triu(row'*ones(1,sz)) +triu(ones(sz,1)*(col.*(col-1))/2);
        ind=triu(ind+skipentries)
        skip=dat.blockind(i);
        dat.mat1=blkdiag(dat.mat1,dat.pack1(ind+ind'-diag(diag(ind))));
        skipentries=skipentries+(sz*(sz+1))/2;
        
      end
      clear ind
      dat=rmfield(dat,'pack1');
    end
  case 'BDFULLV0'
    dat.pack1=fread(fid,dat.pval1,'real*8');
    if(unpack)
      dat.mat1=[];
      skip=0;
      skipentries=0;
      for i=1:dat.nblocks
        
        sz=dat.blockind(i)-skip;
        dat.mat1=blkdiag(dat.mat1,reshape(dat.pack1(skipentries+1:skipentries+sz^2),sz,sz));
        
        %reset indices
        skip=dat.blockind(i);
        skipentries=skipentries+sz^2;
        
      end
      dat=rmfield(dat,'pack1');

    end
 case 'FULLSQV0'

  dat.pack1=fread(fid,dat.pval1^2,'real*8');
  if(unpack)
    dat.mat1=reshape(dat.pack1,dat.nval1,dat.nval1);
      dat=rmfield(dat,'pack1');
  end
  
end 

fclose(fid);

