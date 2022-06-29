function [X,map,out3]=bmpread(filename);
%BMPREAD Read a BMP (Microsoft Windows Bitmap) file from disk.
%	[X,MAP]=BMPREAD('filename') reads the file 'filename' and
%	returns the indexed image X and associated colormap
%	MAP. If no extension is given for the filename, the
%	extension '.bmp' is assumed.
%	[R,G,B]=BMPREAD('filename') reads the 24-bit BMP file
%	from the file 'filename'.
%	BPP=BMPREAD('filename') returns the number of bits per
%	pixel in the BMP file.
%
%	BMPREAD does not read 1-bit or compressed BMP files.
%
%	See also: BMPWRITE, GIFREAD, HDFREAD, PCXREAD, TIFFREAD,
%	          XWDREAD.

%	Mounil Patel 3/10/94
%	Revised Steve Eddins February 9, 1995
%	Copyright (c) 1994 by The MathWorks, Inc.
%	$Revision: 1.13 $  $Date: 1995/02/09 15:17:43 $


if (nargin~=1)
	error('Requires a filename as an argument.');
end;

if (isstr(filename)~=1)
	error('Requires a string filename as an argument.');
end;

if (isempty(findstr(filename,'.'))==1)
	filename=[filename,'.bmp'];
end;

fid=fopen(filename,'rb','l');
if (fid==-1)
	error(['Error opening ',filename,' for input.']);
end;

bfType=fread(fid,2,'uchar');
if (bfType~=[66;77])
	fclose(fid);
	error('Not a BMP file.');
end;

bfSize=fread(fid,1,'uint');
bfReserved1=fread(fid,1,'ushort');
bfReserved2=fread(fid,1,'ushort');
bfOffBytes=fread(fid,1,'uint');

biSize=fread(fid,1,'uint');
if (biSize~=40)
	fclose(fid);
	error('Not a MS Windows Device Independent Bitmap BMP file.');
end;

biWidth=fread(fid,1,'uint');
biHeight=fread(fid,1,'uint');
biPlanes=fread(fid,1,'ushort');
biBitCount=fread(fid,1,'ushort');
if (biBitCount == 1)
  error('BMPREAD does not read 1-bit BMP files.');
end
biCompression=fread(fid,1,'uint');
if biCompression~=0 then
	error('Can''t load compressed format bitmaps.');
end;
biSizeImage=fread(fid,1,'uint');
biXPels=fread(fid,1,'uint');
biYPels=fread(fid,1,'uint');
biClrUsed=fread(fid,1,'uint');
biClrImportant=fread(fid,1,'uint');

if (nargout <= 1)
  X = biBitCount;
  map = [];
  fclose(fid);
  return;
end

if ((nargout == 2) & (biBitCount == 24))
  error('Three output arguments required for 24-bit BMP file.');
end

if ((nargout == 3) & (biBitCount < 24))
  error('Only two arguments needed for 4- and 8-bit BMP files.');
end

if (biBitCount == 24)
  if (rem(biWidth,4)~=0)
    XSize=biWidth+(4-rem(biWidth,4));
  else
    XSize=biWidth;
  end
  YSize=biHeight;
  Size=XSize*YSize;

  fseek(fid, bfOffBytes, 'bof');
  X=fread(fid,Size*3,'uchar');
  out3 = rot90(reshape(X(1:3:length(X)), XSize, YSize)) * (1/255);
  map = rot90(reshape(X(2:3:length(X)), XSize, YSize)) * (1/255);
  X = rot90(reshape(X(3:3:length(X)), XSize, YSize)) * (1/255);
  if (rem(biWidth,4)~=0)
    X=X(:,1:biWidth);
    map = map(:,1:biWidth);
    out3 = out3(:,1:biWidth);
  end;
  
  fclose(fid);
  return;
end  

if (biClrUsed==0)
	nColors=2.^biBitCount;
else
	nColors=biClrUsed;
end;

% load color map now
if (biClrUsed>256)
	map=[];	% 24-bit images have no colormap
else
	map=fread(fid,4*nColors,'uchar');
	map=reshape(map,4,nColors);
	map=map';
	map=map(:,1:3);
	map=fliplr(map);
end;

map=map./255;

% read in 8-bit image data
if (biBitCount==8)
	if (rem(biWidth,4)~=0)
		XSize=biWidth+(4-rem(biWidth,4));
	else
		XSize=biWidth;
	end
	YSize=biHeight;
	Size=XSize*YSize;

	fseek(fid, bfOffBytes, 'bof');
	X=fread(fid,Size,'uchar');
	X=reshape(X,XSize,YSize);
	X=rot90(X);
	X=X+1;
	if (rem(biWidth,4)~=0)
		X=X(:,1:biWidth);
	end;
end;

if (biBitCount==4)
	XSize=ceil(biWidth/2);
	if (rem(XSize,4)~=0)
		XSize=XSize+rem(XSize,4);
	end;
	YSize=biHeight;
	Size=XSize*YSize;
	fseek(fid, bfOffBytes, 'bof');
	X=fread(fid,Size,'uchar');
	X=reshape(X,XSize,YSize);
	X=X';
	loX=X;

	index=loX>127;
	loX(index)=loX(index)-128;
	index=loX>63;
	loX(index)=loX(index)-64;
	index=loX>31;
	loX(index)=loX(index)-32;
	index=loX>15;
	loX(index)=loX(index)-16;

	X=X-loX;
	X=X./16;
	[m,n]=size(X);
	X(:,1:2:(n*2))=X;
	X(:,2:2:(n*2))=loX;
	X=flipud(X);
	X=X+1;
end;

cmax=max(max(X));
map=map(1:cmax,:);

fclose(fid);


