const
 NumJoints = 6;
 NumScrews = 1;

// Global Variables
var
  irobot, iB5, iB6,iscrew,ihead: integer;
  d1, d2, d3: double;
  
  R01, R12, R23,Xw, Rt,mLt,Rlt: Matrix;
  R03: Matrix;

  JointVoltages: array[0..NumJoints - 1] of double;
  JointPos: array[0..NumJoints - 1] of double;
  ScrewH: array[0..NumScrews - 1] of double;

  B6Rot, B6RotCalc, R03Calc: matrix;

function DHMat(a, alpha, d, theta: double): Matrix;
var ct, st, ca, sa: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);
  ca := cos(alpha);
  sa := sin(alpha);
  
  R := Meye(4);
  MSetV(R, 0, 0, ct); MSetV(R, 0, 1,-st * ca);  MSetV(R, 0, 2, st * sa);  MSetV(R, 0, 3, a * ct);
  MSetV(R, 1, 0, st); MSetV(R, 1, 1, ct * ca);  MSetV(R, 1, 2,-ct * sa);  MSetV(R, 1, 3, a * st);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, sa     );  MSetV(R, 2, 2, ca     );  MSetV(R, 2, 3, d     );
  result := R;
end;


function RotZMat(theta: double): Matrix;
var ct, st: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);

  R := Meye(3);
  MSetV(R, 0, 0, ct); MSetV(R, 0, 1,-st);  MSetV(R, 0, 2, 0);
  MSetV(R, 1, 0, st); MSetV(R, 1, 1, ct);  MSetV(R, 1, 2, 0);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, 0 );  MSetV(R, 2, 2, 1);

  result := R;
end;


function RotXMat(theta: double): Matrix;
var ct, st: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);

  R := Meye(3);
  MSetV(R, 0, 0, 1 ); MSetV(R, 0, 1, 0 );  MSetV(R, 0, 2, 0  );
  MSetV(R, 1, 0, 0 ); MSetV(R, 1, 1, ct);  MSetV(R, 1, 2, -st);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, st);  MSetV(R, 2, 2, ct );

  result := R;
end;

 function R36(teta4, teta5, teta6: double): matrix;
  var s4, s5, s6, c4, c5, c6: double;
  begin

    s4:=sin(teta4);
    s5:=sin(teta5);
    s6:=sin(teta6);
    c4:=cos(teta4);
    c5:=cos(teta5);
    c6:=cos(teta6);
    result:=Mzeros(3,3);
    Msetv(result, 0, 0, c5);
    Msetv(result, 0, 1, -s5*c6);
    Msetv(result, 0, 2, s5*s6);
    Msetv(result, 1, 0, c4*s5);
    Msetv(result, 1, 1, c4*c6*c5-(s4*s6));
    Msetv(result, 1, 2, -s4*c6-(c4*c5*s6));
    Msetv(result, 2, 0, s5*s4);
    Msetv(result, 2, 1, c4*s6+c5*c6*s4);
    Msetv(result, 2, 2, c4*c6-s4*c5*s6);


  end;

    function R03T(teta1, teta2, teta3: double): matrix;
  var s1, s2, s3, c1, c2, c3, teta4, teta5, teta6: double;
  begin

    s1:=sin(teta1);
    s2:=sin(teta2);
    s3:=sin(teta3);
    c1:=cos(teta1);
    c2:=cos(teta2);
    c3:=cos(teta3);
    result:=Mzeros(3,3);
    Msetv(result, 0, 0, c1*c2*c3-c1*s2*s3);
    Msetv(result, 0, 1, c2*c3*s1-s1*s2*s3);
    Msetv(result, 0, 2, -s2*c3-c2*s3);
    Msetv(result, 1, 0, -s1);
    Msetv(result, 1, 1, c1);
    Msetv(result, 1, 2, 0);
    Msetv(result, 2, 0, c1*c2*s3+c1*s2*c3);
    Msetv(result, 2, 1, c2*s1*s3+c3*s1*s2);
    Msetv(result, 2, 2, -s2*s3+c2*c3);

    teta4:=atan2(c1*c2*s3+c1*s2*c3, -s1);
    teta6:=atan2(-s2*c3-c2*s3, -(c2*c3*s1-s1*s2*s3));
    teta5:=atan2(((-s2*c3-c2*s3)/(sin(teta6))), c1*c2*c3-c1*s2*s3);

    SetAxisPosRef(irobot,3,teta4);
    SetAxisPosRef(irobot,4,teta5);
    SetAxisPosRef(irobot,5,teta6);

    SetRCValue(27, 5, format('%.3g',[deg(teta4)]));
    SetRCValue(28, 5, format('%.3g',[deg(teta5)]));
    SetRCValue(29, 5, format('%.3g',[deg(teta6)]));

  end;

function mxtrans(xt, yt, zt: double): matrix;
  begin

  result:=Mzeros(3,1);
  Msetv(result, 0, 0, xt);
  Msetv(result, 1, 0, yt);
  Msetv(result, 2, 0, zt);
  end;

  function mLtrans(Lt: double): matrix;
  begin

    result:=Mzeros(3,1);
    Msetv(result, 0, 0, Lt);
    Msetv(result, 1, 0, 0);
    Msetv(result, 2, 0, 0);

  end;


function IK(Xtt,Rt:matrix;Lt:double):matrix;
var xt,yt,zt:double;
    updown, dt, h, dh, alfa3, alfa3num, alfa3den, alfa3aux, alfah, alfahaux, tetat: double;
    teta:array[1..3] of double;
begin

    mLt:=mLtrans(Lt);
    Rlt:= Mmult(Rt, mLt);
    Xw:= Msub(Xtt, Rlt);
    MatrixToRangeF(11, 7, Xtt, '%.3f');
    MatrixToRangeF(11, 6, Xw, '%.3f');

    xt:= Mgetv(Xw,0,0);
    yt:= Mgetv(Xw,1,0);
    zt:= Mgetv(Xw,2,0);


    if ((xt>1) or (xt<-1) or (yt>1) or (yt<-1) or (zt>1) or (zt<-1)) then begin
      SetRCValue(9, 4, 'Erro nos valores');
      exit;
    end else begin

    updown:=GetRCValue(4, 4);

    if ((updown<>1) and (updown<>0)) then begin
    SetRCValue(4, 5, 'Erro up/down');
    exit;
    end;

    if(updown=1) then begin

    dt:=sqrt(sqr(xt)+sqr(yt));
    teta[1]:=atan2(yt, xt);
    h:=zt-d1;
    dh:=sqrt(sqr(dt)+sqr(h));
    alfa3num:=(sqr(d3)-sqr(d2)-sqr(dh));
    alfa3den:=(-2*d2*dh);
    if (alfa3den=0) then begin
      SetRCValue(9, 4, 'Erro den alfa3=0');
      exit;
      end;
    alfa3aux:=(alfa3num/alfa3den);
    if ((alfa3aux>1) or (alfa3aux<-1)) then begin
      SetRCValue(9, 4, 'Erro arcos alfa3');
      exit;
    end;
    alfa3:=arccos(alfa3aux);
    tetat:=atan2(h, dt);
    teta[2]:=-(alfa3+tetat);
    alfahaux:=((sqr(dh)-sqr(d2)-sqr(d3))/(-2*d2*d3));
    if ((alfahaux>1) or (alfahaux<-1)) then begin
      SetRCValue(9, 4, 'Erro arcos alfah');
      exit;
    end;
    alfah:=arccos(alfahaux);
    teta[3]:=pi-alfah;



    SetRCValue(24, 5, format('%.3g',[deg(teta[1])]));
    SetRCValue(25, 5, format('%.3g',[deg(teta[2])]));
    SetRCValue(26, 5, format('%.3g',[deg(teta[3])]));

    SetAxisPosRef(irobot,0,teta[1]);
    SetAxisPosRef(irobot,1,teta[2]);
    SetAxisPosRef(irobot,2,teta[3]);

      R03Calc:=R03T(teta[1], teta[2], teta[3]);
      MatrixToRangeF(21, 4, R03Calc, '%.3f');

    end;

    if (updown=0) then begin

     dt:=sqrt(sqr(xt)+sqr(yt));
    teta[1]:=atan2(yt, xt);
    h:=zt-d1;
    dh:=sqrt(sqr(dt)+sqr(h));
    alfa3num:=(sqr(d3)-sqr(d2)-sqr(dh));
    alfa3den:=(-2*d2*dh);
    if (alfa3den=0) then begin
      SetRCValue(9, 4, 'Erro den alfa3=0');
      exit;
      end;
    alfa3aux:=(alfa3num/alfa3den);
    if ((alfa3aux>1) or (alfa3aux<-1)) then begin
      SetRCValue(9, 4, 'Erro arcos alfa3');
      exit;
    end;
    alfa3:=arccos(alfa3aux);
    tetat:=atan2(h, dt);
    teta[2]:=(alfa3-tetat);
    alfahaux:=((sqr(dh)-sqr(d2)-sqr(d3))/(-2*d2*d3));
    if ((alfahaux>1) or (alfahaux<-1)) then begin
      SetRCValue(9, 4, 'Erro arcos alfah');
      exit;
    end;
    alfah:=arccos(alfahaux);
    teta[3]:=alfah-pi;

    SetRCValue(24, 5, format('%.3g',[deg(teta[1])]));
    SetRCValue(25, 5, format('%.3g',[deg(teta[2])]));
    SetRCValue(26, 5, format('%.3g',[deg(teta[3])]));

    SetAxisPosRef(irobot,0,teta[1]);
    SetAxisPosRef(irobot,1,teta[2]);
    SetAxisPosRef(irobot,2,teta[3]);

      R03Calc:=R03T(teta[1], teta[2], teta[3]);
      MatrixToRangeF(21, 4, R03Calc, '%.3f');

    end;




  end;








end;


procedure UpdateScrew(index: integer);
var w, h: double;
begin
  w := GetAxisSpeed(index, 1); // second axis has the head rotation
  SetRCValue(3, 4, format('%.2g',[w]));
  if abs(w) > 0.1 then begin
    ScrewH[index - 1] := ScrewH[index - 1] + w/1000;
    SetAxisPosRef(index, 0, ScrewH[index - 1]);
  end;

end;


procedure Control;
var i: integer;
    B5Pos, B6Pos, Rt, Xw, Rlt, mxt: Matrix;
    headPos, headRot: Matrix;
    //xw, yw, zw: double;
    x, y, z, d, updown, dt, h, dh, alfa3, alfa3num, alfa3den, alfa3aux, alfah, alfahaux, tetat, Lt: double;
    RTool, XTool, XWrist, LTool: matrix;
    teta:array[1..3] of double;
begin
  UpdateScrew(1);

  B5Pos := GetSolidPosMat(iRobot, iB5);
  MatrixToRange(11, 2, B5Pos);

  B6Pos := GetSolidPosMat(iRobot, iB6);
  MatrixToRange(16, 2, B6Pos);

  B6rot := GetSolidRotMat(iRobot, iB6);
  MatrixToRangeF(16, 4, B6Rot, '%.3f');



  // Read joint positions
  for i := 0 to NumJoints -1 do begin
    JointPos[i] := GetAxisPos(irobot, i);
  end;

  // and show
  for i := 0 to NumJoints -1 do begin
    SetRCValue(3 + i, 2, format('%.3g',[Deg(JointPos[i])]));
  end;


   d:=d2*cos(-JointPos[1])+d3*cos(-JointPos[1] -JointPos[2]);
  x:=d*cos(JointPos[0]);
  y:=d*sin(JointPos[0]);
  z:=d1+ d2*sin(-JointPos[1])+d3*sin(-JointPos[1]-JointPos[2]);
  SetRCValue(11, 3, format('%.3g',[x]));
  SetRCValue(12, 3, format('%.3g',[y]));
  SetRCValue(13, 3, format('%.3g',[z]));
  SetRCValue(3, 4, 'Up');
  B6Rot:=getSolidRotMat(iRobot, iB6);
  MatrixToRangeF(16, 4, B6Rot, '%.3f');
  B6RotCalc:=R36(JointPos[3], JointPos[4], JointPos[5]);
  MatrixToRangeF(16, 8, B6RotCalc, '%.3f');
  SetRCValue(24, 4, 'teta1:');
  SetRCValue(25, 4, 'teta2:');
  SetRCValue(26, 4, 'teta3:');
  SetRCValue(27, 4, 'teta4:');
  SetRCValue(28, 4, 'teta5:');
  SetRCValue(29, 4, 'teta6:');
  SetRCValue(5, 6, 'Rt');
  SetRCValue( 10, 6, 'Xw');

  // control equations
  Rt:=RangeToMatrix(6, 6, 3, 3);
  Lt:=0.035;



  if RCButtonPressed(10,4) then begin

    //xt:=GetRCValue(11,4);
    //yt:=GetRCValue(12,4);
   // zt:=GetRCValue(13,4);
    headPos:=getSolidPosMat(iscrew,ihead);
    headRot:=getSolidRotMat(iscrew,ihead);


    IK(headPos,headRot,Lt);

  end;
  // ...

end;

procedure Initialize;
var i: integer;
begin
  irobot := 0;
  iscrew:= 1;
  ihead:= getSolidIndex(iscrew,'screw_head');
  iB5 := GetSolidIndex(irobot, 'B5');
  iB6 := GetSolidIndex(irobot, 'B6');
  SetRCValue(2, 1, 'Joint');
  SetRCValue(2, 2, 'Pos (deg)');
  for i := 0 to NumJoints -1 do begin
    SetRCValue(3 + i, 1, format('%d',[i]));
  end;

  for i := 0 to NumScrews -1 do begin
    Screwh[i] := GetAxisPos(1 + i, 0);
  end;

  d1 := 0.55;
  d2 := 0.4;
  d3 := 0.37;
end;
