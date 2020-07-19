    // Motion is a vectorField of all moving boundary points
    vectorField motion(ms.movingPoints().size(), vector::zero);

    vectorField oldPoints = ms.movingPoints();

    scalar oldTime = time_.value() - time_.deltaT().value();
    scalar curTime = time_.value();

    scalar alphaOld   = 0.0;
    scalar alphaCur   = 0.0;

//     scalar rotationAmplitude_(0.707);
//     scalar rotationFrequency_(0.0);
//     vector translationFrequency_(0.0,0.25,0.0);
//     vector translationAmplitude_(0.0,5.0,0.0);
//     vector initialRotationOrigin_(0.0,0.0,0.0);

    scalar pi = 3.141592653589764;
//    scalar phi = pi/2;
    scalar aoaOld = 0.0;
    scalar aoaCur = 0.0;
    scalar ydotOld = 0.0;
    scalar ydotCur = 0.0;
    scalar xdotOld = 0.0;
    scalar xdotCur = 0.0;
    scalar Tp = 1.0/rotationFrequency_;
    scalar THalf = Tp/2.0;

    scalar magU = 1.1;
    scalar chord = 1.0;

    xdotOld = -2.0*pi*translationFrequency_[0]/Foam::tan(beta_)*Foam::sin(2.0*pi*translationFrequency_[0]*oldTime);
    xdotCur = -2.0*pi*translationFrequency_[0]/Foam::tan(beta_)*Foam::sin(2.0*pi*translationFrequency_[0]*curTime);

    ydotOld = -2.0*pi*translationFrequency_[1]*Foam::sin(2.0*pi*translationFrequency_[1]*oldTime);
    ydotCur = -2.0*pi*translationFrequency_[1]*Foam::sin(2.0*pi*translationFrequency_[1]*curTime);

//    aoaOld = rotationAmplitude_*Foam::sin(2.0*pi*rotationFrequency_*oldTime);
//    aoaCur = rotationAmplitude_*Foam::sin(2.0*pi*rotationFrequency_*curTime);

    if (curTime < THalf)
	{
	    aoaCur = rotationAmplitude_*(0.5-0.5*Foam::cos(4.0*pi*rotationFrequency_*curTime));
	}
    else if ((curTime >= THalf) && (curTime < Tp))
	{
	    aoaCur = 0.0;
	}
    else if ((curTime >= Tp) && (curTime < 1.5*Tp))
	{
	    aoaCur = rotationAmplitude_*(0.5-0.5*Foam::cos(4.0*pi*rotationFrequency_*curTime)); 
	}
    else if ((curTime >= 1.5*Tp) && (curTime < 2.0*Tp))
	{
	    aoaCur = 0.0;
	}

    if (oldTime < THalf)
	{
	    aoaOld = rotationAmplitude_*(0.5-0.5*Foam::cos(4.0*pi*rotationFrequency_*oldTime));
	}
    else if ((oldTime >= THalf) && (oldTime < Tp))
	{
	    aoaOld = 0.0;
	}
    else if ((oldTime >= Tp) && (oldTime < 1.5*Tp))
	{
	    aoaOld = rotationAmplitude_*(0.5-0.5*Foam::cos(4.0*pi*rotationFrequency_*oldTime)); 
	}
    else if ((oldTime >= 1.5*Tp) && (oldTime < 2.0*Tp))
	{
	    aoaOld = 0.0;
	}
   

//    alphaOld = (Foam::atan(ydotOld/(magU + xdotOld)) - aoaOld);
//    alphaCur = (Foam::atan(ydotCur/(magU + xdotCur)) - aoaCur);

    alphaOld = - (Foam::atan2(ydotOld,(magU - xdotOld)) + aoaOld);
    alphaCur = - (Foam::atan2(ydotCur,(magU - xdotCur)) + aoaCur);

//    alphaOld = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*oldTime);
//    alphaCur = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*curTime);


    vector translationVector
    (
        translationAmplitude_[0]/Foam::tan(beta_)*
        (
            (Foam::cos(2*pi*translationFrequency_[0]*curTime) + chord)
          - (Foam::cos(2*pi*translationFrequency_[0]*oldTime) + chord)
        ),
        translationAmplitude_[1]*
        (
            (Foam::cos(2*pi*translationFrequency_[1]*curTime) - chord)
          - (Foam::cos(2*pi*translationFrequency_[1]*oldTime) - chord)
        ),
        0
    );

    tensor RzOld
    (
        Foam::cos(alphaOld), -Foam::sin(alphaOld), 0,
        Foam::sin(alphaOld),  Foam::cos(alphaOld), 0,
        0, 0, 1
    );

    tensor RzCur
    (
        Foam::cos(alphaCur), -Foam::sin(alphaCur), 0,
        Foam::sin(alphaCur),  Foam::cos(alphaCur), 0,
        0, 0, 1
    );

    vectorField rotationField
    (
        (RzCur - RzOld)
      & (statPoints_ - initialRotationOrigin_)
    );

    motion = translationVector + rotationField;
