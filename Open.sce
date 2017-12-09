/*
Generation of a three-dimensional GMSH .geo file
from a BlueKenue two-dimensional surface depicting 
a DEM topography. 

Blocks: 
1) functions definition (FDD)
2) topography face construction (TFC) 
3) bottom base face construction (BBC)
4) stitching of up&down boundaries (SBS)
5) create CVS file to export (CVS)
*/


//*********************FDD*********************

//Define Cantor Pairing Function
clear
clc

function [c] = cantorP(a,b)
    c=(((a+b)*(a+b+1))/2)+b;
endfunction

//Read point data and elements connectivities
[points]=csvRead('./EMSH/points.csv');
[elements]=csvRead('./EMSH/elements.csv');

//*********************************************

//*********************TFC*********************

//Determine points and elements from Bluesurface mesh
nPoints=size(points,1);
nElements=size(elements,1);
nDimensions=size(points,2);
nIJ=[1 2 3 1];

//Generate list of points GMSH_format
for i = 1:nPoints
    temp = strcat(['Point(',string(i),') = {',string(points(i,1)),',',string(points(i,2)),',',string(points(i,3)),',ref};']);
    gmshPoints(i)= temp;
end
clear temp

//Build unique Id for lines
for r=1:nElements
    for x= 1:size(nIJ,2)-1
        if elements(r,nIJ(x))>elements(r,nIJ(x+1)) then
            elements(r,x+3)=cantorP(elements(r,nIJ(x)),elements(r,nIJ(x+1)));
        else
            elements(r,x+3)=cantorP(elements(r,nIJ(x+1)),elements(r,nIJ(x)));
        end
    end
end

//Generate list of lines GMSH_format
A=[];
c=0;
for l = 1:nElements
    if find(A==elements(l,4))==[] then
        temp = strcat(['Line(',string(elements(l,4)),') = {',string(elements(l,1)),',',string(elements(l,2)),'};']);
        c=c+1;
        gmshLines(c)=temp;
        A(l,1)=elements(l,4);
    else
        A(l,1)=-elements(l,4);
    end

    if find(A==elements(l,5))==[] then
        temp = strcat(['Line(',string(elements(l,5)),') = {',string(elements(l,2)),',',string(elements(l,3)),'};']);
        c=c+1;
        A(l,2)=elements(l,5);
        gmshLines(c)=temp;
    else
        A(l,2)=-elements(l,5);
    end

    if find(A==elements(l,6))==[] then
        temp = strcat(['Line(',string(elements(l,6)),') = {',string(elements(l,3)),',',string(elements(l,1)),'};']);
        c=c+1;
        A(l,3)=elements(l,6);
        gmshLines(c)=temp;
    else
        A(l,3)=-elements(l,6);
    end
end

//Build Line Loops & plane Surface GMSH_format
llId=max(elements(:,4:6));

for r = 1:nElements
    temp = strcat(['Line Loop(',string(llId+r),') = {',string(A(r,1)),',',string(A(r,2)),',',string(A(r,3)),'};']);
    gmshLoops(r)= temp;
    temp2 = strcat(['Plane Surface(',string(llId+nElements+r),') = {',string(llId+r),'};']);
    gmshSurface(r)= temp2;
    S(r,1)=(llId+nElements+r);
end

//*********************************************


//*********************BBC*********************

//Extract perimeter and identify blines
n=0;
for r=1:size(A,1)
    for c=1:nDimensions
        if find(-A(r,c)==A)==[] then
            n=n+1;
            //Lines
            bL(n,1)=A(r,c);
            bL(n,2)=elements(r,nIJ(c));
            bL(n,3)=elements(r,nIJ(c+1));
            //Coordinates
            bL(n,4)=points(bL(n,2),1);
            bL(n,5)=points(bL(n,2),2);
            bL(n,6)=0;
        end
    end
end

//Build bottom points
for r = 1:size(bL,1)
    temp = strcat(['Point(',string(nPoints+bL(r,2)),') = {',string(bL(r,4)),',',string(bL(r,5)),',',string(bL(r,6)),',ref};']);
    gmshBPoints(r)= temp;
end

//Build bottom lines GMSH_format
c=1;
lbId = llId + (2*nElements)+2;

for r = lbId+1:lbId+size(bL,1)
    temp = strcat(['Line(',string(lbId+bL(c,1)),') = {',string(nPoints+bL(c,2)),',',string(nPoints+bL(c,3)),'};']);
    c=c+1;
    gmshBLines(c)=temp;
end

//Build bottom loop

B(1,1)=bL(1,1);
B(1,2)=bL(1,2);
B(1,3)=bL(1,3);
for r = 2:size(bL,1)
    B(r,2)=B(r-1,3);
    m = find(bL(:,2)==B(r,2));
    B(r,1)=bL(m,1);
    B(r,3)=bL(m,3);
end

//Build bottom Loop GMSH_format
B(:,4)=B(:,1);
B(:,1)=B(:,1)+lbId;
lcId = max(B(:,1));

hemp = strcat(['Line Loop(',string(lcId+1),') = {']);
temp =[];
for r=1:size(B,1)
    if r == size(B,1) then
        temp = strcat([temp,string(B(r,1))]);
    else
        temp = strcat([temp,string(B(r,1)),',']);
    end
end
gmshBLoops=strcat([hemp,temp,'};']);

//Build bottom Surface GMSH_format
temp = strcat(['Plane Surface(',string(lcId+2),') = {-',string(lcId+1),'};']);
gmshBSurface = temp;

//*********************************************

//*********************SBS*********************

ldId = max(B(:,1));+1;
for r = 1:size(B,1)
    temp = strcat(['Line(',string(ldId+r),') = {',string(B(r,2)),',',string(B(r,2)+nPoints),'};']);
    gmshCLines(r)=temp;
    C(r)=ldId+r;
end
C(size(B,1)+1)=C(1)
clear temp

leId = ldId + size(B,1);
for r = 1:size(B,1)
    temp = strcat(['Line Loop(',string(leId+r),') = {',string(C(r)),',',string(B(r,1)),',-',string(C(r+1)),',-',string(B(r,4)),'};']);
    gmshCLoops(r)=temp;
    hemp = strcat(['Plane Surface(',string(leId+size(B,1)+r),') = {',string(leId+r),'};']);
    gmshCSurface(r)=hemp;
    S2(r,1)=leId+size(B,1)+r;
end


temp=[];
for r = 1:size(S,1)+size(S2,1)+1
    if r <= size(S,1) then
        if r == size(S,1) then
            temp = strcat([temp,',',string(S(r,1))]);
        elseif r == 1 then
            temp = strcat([string(S(r,1))]);
        else
            temp = strcat([temp,',',string(S(r,1))]);
        end
    elseif r == size(S,1)+size(S2,1)+1 then
        temp = strcat([temp,',',string(lcId+2)]);
    else   
        if r == size(S,1)+size(S2,1) then
            temp = strcat([temp,',',string(S2(r-size(S,1),1))]);
        elseif r == size(S2,1) then
            temp = strcat([string(S2(r-size(S,1),1))]);
        else
            temp = strcat([temp,',',string(S2(r-size(S,1),1))]);
        end
    end
gmshSurfLoop=strcat(['Surface Loop(',string('1'),') ={',temp,'};']);
end
clear temp

//*********************************************

GMSHV = strcat(['Volume(1) = {1};']);

//Write .geo GMSH file

GMSHA = ['ref=100;' ; gmshPoints ; gmshLines ; gmshLoops ; gmshSurface];
GMSHB = ['ref=500;' ; gmshBPoints ; gmshBLines ; gmshBLoops ; gmshBSurface];
GMSHC = [gmshCLines ; gmshCLoops ; gmshCSurface ; gmshSurfLoop];
GMSH = [GMSHA ; GMSHB ; GMSHC ; GMSHV];
csvWrite(GMSH,'./EMSH/gmsh.geo');
