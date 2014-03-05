function alphaRadiusVector = getAlphaRadius( dataStr )
% this function returns the alpha shape radii for each time step
% of the current data set; the radii are precomputes as an optimal value in
% scifer using the CGAL lib
if strcmp( dataStr, '120830_raw' )
  alphaRadiusVector = [
    47013.4;48351.7;47784.1;47225.3;46676.2;46137.6;45610.4;45095;44592.2;44102.3;...
    43625.9;43163.2;42714.6;42280.3;41860.5;41455.3;41064.8;40689.2;39898.5;33857;...
    32954.6;40634.6;49259.3;34484.9;35180.6;35922.8;36715.6;37563.8;38472.9;39448.7;...
    40498.1;41628.8;42849.5;44170.3;45602.6;38182.2;39478.6;40891.4;42435.2;44127.1;...
    45987;48038.6;50309.9;36507.1;37730;39040.4;40447.1;41960;43590.4;45351;...
    47256.7;49323.7;51572.3;4665.65;4636.53;4602.48;2609.75;2473.4;2344.32;2466.31;...
    2644.55;2832.16;3029.4;3236.51;3453.74;3681.32;3919.51;4168.53;4428.64;4700.05;...
    4983.03;5277.78;5584.56;5903.59;6235.09;6579.3;6936.43;7306.72;7690.37;8087.62;...
    8425.7;8510.74;8595.95;8681.54;8835.84;9309.26;9801.32;10312.5;10843.3;11394.2;...
    11965.7;11950.8;11785.3;11624.1;11466.9;11313.6;11164.2;11018.3;10875.9;10736.9;...
    10601.1;10468.4;10338.8;10270.8;10298.1;9966.94;9880.65;9868.57;9854.54;9838.58;...
    4002.44;3993.12;3984.22;3975.76;3967.71;3960.06;3952.8;3945.92;3939.4;3933.25;...
    3927.44;3921.98;3916.83;3912.01;1324.13;1332.9;1341.74;1350.66;1359.65;1368.71;...
    1377.85;1387.06;1396.35;1405.71;1447.6;1495.5;1546.8;1601.5;1659.57;1720.99;...
    1785.74;1853.78;1925.1;1999.67;2077.47;2158.47;2242.64;2329.96;2420.4;2513.94;...
    2610.55;2710.2;2812.86;2918.51;3027.12;3138.67;3253.12;3370.45;3490.63;3613.63;...
    3739.43;3868;3999.31;4133.35;4270.07;4409.44;4551.46;4696.08;4843.29;4993.05;...
    5145.35;5300.15;5457.43;5617.17;5779.33;5881.26;5943.39;6004.35;6064.14;6122.75;...
    6180.16;6236.39;6291.42;6345.24;6397.87;6449.3;6499.52;6548.54;6596.37;6643;...
    6688.44;6732.71;6775.79;6817.7;6858.44;6898.03;6936.47;6973.76;7009.93;7044.97;...
    7078.9;7111.73;7143.46;7174.11;7203.69;7232.21;7259.68;7286.11;7311.52;7335.91;...
    7359.29;7381.69;7403.12;7423.57;7443.07;7461.63;7479.25;7495.97;7511.78;7526.7;...
    7540.74;7553.92;7566.24;7577.72;7588.38;7598.22;7607.26;7615.52;7622.99;7629.7;...
    7635.66;7640.88;7645.38;7649.17;7652.25;7654.65;7634.25;7565.16;7497.21;7430.37;...
    7364.63;7299.93;7236.29;7173.66;7112.02;7051.36;6991.65;6932.87;6875;6818.03;...
    6761.93;6706.69;5574.42;5161.25;5040.97;5148.78;5266.06;5363.4;5440.9;5498.82;...
    5537.57;5557.67;5559.77;5544.59;5512.9;5465.56;5403.42;5327.39;5238.41;5137.44;...
    5025.44;4903.39;4772.26;4633.05;4486.75;4334.34;4176.8;4015.11;3850.25;3683.19;...
    3514.88;3346.27;3178.3;3011.91;2882.41;2766.84;2666.6;2581.8;2512.53;2458.93;...
    2421.12;2264.88;2140.39;2050.04;1987.27;1946.47;1922.8;1912.2;1911.21;1916.95 ];
elseif strcmp( dataStr, '121204_raw_2014' )
  alphaRadiusVector = [
    49357.7;52505.9;52123.3;51760.6;51416.9;51091.6;50071.3;48987.6;47916.2;46856.1;...
    45806.3;44765.6;43733.7;42709.5;38806.6;37933.6;37058.3;36181.3;35302.9;34423.5;...
    33544;32664.5;31785.6;30907.9;30031.9;29158;28287.1;27505.7;27915.8;28330.8;...
    28750.9;29175.9;29606.1;30041.3;30481.6;30927.2;38820.5;34644.5;32295.2;32761.9;...
    33233.9;33711.3;34194.2;34682.6;35176.6;35676.1;36181.2;36692;37208.6;37730.8;...
    38258.9;38792.8;39332.6;39878.3;40430;40987.7;41551.4;42121.3;42697.3;43279.5;...
    43867.9;44462.6;45063.7;45671.1;46285;46905.3;47532.1;48165.5;48805.5;49452.1;...
    50105.5;50765.6;51432.4;52106.2;52786.8;53474.4;54168.9;54870.5;46832.4;46917.3;...
    47003.3;47090.4;47178.5;47267.9;47358.6;47450.6;47544.1;47639;47735.5;47833.7;...
    47933.5;48035.1;48138.5;48243.9;48351.1;48460.5;48572;48685.6;48801.5;48919.7;...
    49040.4;49163.4;49717.8;49740.4;49763.7;49787.7;49812.5;49838.1;49864.5;49891.9;...
    49920.2;49949.6;49980;50011.6;50044.3;50040.1;49878.2;49722.7;49573.7;49431.5;...
    49296.3;49168.3;49047.8;48935;48830.2;51524.3;51638.6;51755.4;51874.9;50711.5;...
    49460.6;48221.7;46995.3;45781.1;44579.7;43390.7;42214.5;41050.9;39900.4;38762.8;...
    37638.5;36527.6;35429.8;34345.9;33275.3;32218.9;31281;30654;30032.3;29416;...
    28805.2;28199.9;27600.3;27006.5;26418.5;25836.4;25260.3;24690.3;24126.5;23569;...
    23017.8;22473;21934.6;21402.9;20877.7;20359.2;19847.4;19342.5;18844.4;18353.3;...
    17869.1;17392;16921.9;16459;16003.2;15554.7;15113.4;14679.5;14252.9;13833.7;...
    13421.9;13017.6;12620.7;12424.2;12358.7;12294.8;12232.5;12171.6;12112;12053.7;...
    11996.5;11940.4;11885.3;11831.1;11777.7;11725;11647.7;11518.6;11389.8;11261.6;...
    10134.6;10027.2;9923.64;9823.75;9696.95;9557.3;5343.31;5502.5;5576.89;5615.15;...
    5653.15;5690.9;5728.43;5745.76;5735.4;5725.51;5716.08;5707.11;5698.6;5690.53;...
    5682.92;5675.74;5669;5662.69;5656.81;5865.94;6339.49;6798.16;7249.23;5633.7;...
    5822.91;6030.61;6243.3;6460.95;6683.58;6911.2;7143.85;7381.51;7624.22;7871.97;...
    8124.78;8382.67;8645.65;8900.82;8738.06;8577.8;8420;8264.62;8111.65;7961.06;...
    7812.83;7668.22;7901.53;8139.84;8383.17;8631.59;8885.13;9143.85;9407.83;9677.16;...
    9951.88;10232.1;10517.9;10809.4;11106.8;11410;11719.3;12034.9;12356.8;12685.2;...
    13020.3;13362.3;13711.4;14067.8;14431.8;14803.5;15183.3;15571.4;43310;16373.7;...
    16788.6;17213.1;17625.5;17592.2;17563.1;17538.4;17518.1;17502.3;17491.1;17276.1;...
    16868.7;16469.6;16079.1;15697.5;10762.6;9769.36;12405.4;16058.5;20703.1;22883.4 ];
elseif strcmp( dataStr, '121211_raw' )
  alphaRadiusVector = [
    18200.7;18510.6;18827.7;19152.3;19484.6;128777;125894;125357;122363;119379;...
    116415;113479;110576;107710;104886;102107;99372.8;96687;94050.2;91463.4;...
    88927.4;86442.4;84008.7;81626.4;79295.4;77015.4;74786.4;72607.9;70479.7;68401.3;...
    66372.3;64392;62460.2;60576.1;58739.2;56949;55204.8;52268.5;48939.1;45780.7;...
    42791.5;39969.3;37312.1;34817.8;32484.1;30308.6;26143.8;24827.7;23886.5;23141.5;...
    21721.1;20444.1;19308.2;18310.5;17448.6;16719.9;16121.6;15651.2;15306.1;15083.5;...
    14980.8;14995.4;15124.5;15365.6;15716;16173;16734;17396.4;18157.5;19014.7;...
    19965.3;21006.8;22136.6;23352.1;24650.6;26029.8;27486.9;29019.4;30625;32301;...
    34044.8;35854.2;37726.6;39659.6;41650.7;43697.6;45797.8;47949.1;50148.9;52395;...
    54685.2;57017.1;54748.6;57134.1;59542.1;61969.8;64415;66875.1;69347.9;71831.1;...
    74322.6;76820;79321.6;81825;84328.5;86830;89327.8;91819.9;94304.8;96780.7;...
    99245.9;101699;104138;106562;108970;111359;113729;116078;118406;120711;...
    122992;125247;37652;67522.5;69778.7;72142.5;74618.7;77211.1;79924;82762.4;...
    85731.1;88836.4;92083.2;95477.1;99025.4;102734;106611;110662;114897;119323;...
    123949;128784;133840;139125;144651;139772;151551;123958;116936;110432;...
    104400;98800.2;93592.8;88747.6;84236;80028.8;76104;72439.7;69017.2;65818.8;...
    62828.6;60032.4;57417.3;54970.7;52682.1;50541.9;48540.8;46671.2;44924.4;43294;...
    41773.9;40358.1;39041.3;37819;36686;35638.8;34673.4;33786.3;32974.7;32235.5;...
    31566;30963.9;30426.9;29953.3;29541.2;29188.8;28895;28658.3;28618.4;29192;...
    29766.7;30342.7;30919.2;31496.1;32073;32649.7;33226.1;33801.2;34375.3;34947.9;...
    35518.6;36087.5;36653.6;37217;37777.3;38334.4;38887.9;39436.7;39981.1;40520.7;...
    41055.1;41584.1;42106.8;42623.5;43133.6;44103.9;45532.1;47035.2;48617.1;50280.2;...
    52027.1;53861.7;55786.4;57805.4;59922.8;62141.7;64467.7;66904.3;69457;72131;...
    74931;77864.9;80937.4;84155.9;87527.1;91059.1;94760.6;98639.4;102705;106967;...
    111437;116127;121048;126213;131637;137333;143322;149617;156238;163205;...
    170540;178271;186416;195005;204070;213640;223755;234453;245768;257752;...
    270453;283924;298229;313422;329583;346788;359090;353701;348468;343386;...
    338447;333646;328980;324442;320028;315733;311552;307483;303520;262059;...
    263865;265685;267519;269367;271229;273106;274996;276902;277852;245874;...
    244009;242171;240360;238576;236818;235085;233377;231694;230035;228400 ];
elseif strcmp( dataStr, '130508_raw' )
  alphaRadiusVector = [
    151466;145698;140116;134717;129494;124444;119560;114838;110274;105863;...
    101601;97484.4;93508.1;89668.7;103208;100555;97956.6;95412.1;92920.6;90481.5;...
    88094.5;85758.4;83472.8;81237.3;79050.9;76913.1;74823.4;72780.9;70785.2;68835.8;...
    66931.7;65072.6;63258;61486.9;59759.1;58073.9;56430.6;51104.6;50956.7;50803.1;...
    50643.8;50478.9;50308.7;50133.3;49952.6;49767;49576.4;49381.1;45399.3;45033.1;...
    44656.7;44270.2;43874.1;43468.5;43053.7;42630;42197.7;41757;41308.2;40430.6;...
    40065.6;39696.9;39324.5;38948.6;38569.1;497545;500518;503184;505553;507622;...
    509396;510882;512088;513011;513663;514049;514175;514048;513674;513060;...
    512219;511153;509871;508383;506700;504825;502769;407400;405196;402909;...
    400551;398130;395652;393131;390571;387980;385369;382746;380117;377492;...
    374879;372286;369721;367191;364704;362269;319675;313952;308308;302745;...
    297263;291865;286553;281327;276190;271141;266182;261313;256536;251850;...
    247257;242755;238348;234033;229810;225681;221645;217702;213850;210091;...
    206424;202848;199362;195967;192661;189444;186315;183273;180317;177447;...
    174662;171960;169341;166803;164346;161969;159670;157449;155303;153233;...
    151237;149313;147461;145679;143967;142323;140745;139233;137785;136400;...
    135077;133815;132612;131467;130379;129347;128369;127444;126572;125751;...
    124979;124256;123580;122951;122366;121825;121328;120871;120456;120080;...
    119742;119442;119178;118949;118754;118593;118463;118365;118297;118259;...
    118248;118265;118309;118378;118471;118589;118729;118892;119075;119279;...
    119503;119746;120006;120284;120579;120889;121214;121554;121907;122274;...
    122652;123043;123444;123856;124278;124709;125148;125596;126051;126513;...
    126981;127456;127935;128419;128908;129401;129897;130396;130897;131401;...
    131906;132412;132919;133427;133934;134441;134948;135453;135957;136459;...
    136959;137457;137951;138443;138932;139417;139897;140374;140846;141314;...
    141777;142234;142686;143132;143573;144007;144435;144856;145271;145679;...
    146080;146474;146860;147239;147610;147973;148328;148675;149014;149344;...
    145583;141597;137719;133945;130272;126697;123218;119832;116535;113325;...
    110200;107157;104193;101308;98497.1;95760;93093.9;90497.5;87968.2;85504.3;...
    83104.4;80766.7;78489.4;76270.7;74109.8;72004.3;69953.5;67955.6;66009.5;64113.7;...
    62267;60468.3;58716.4;57010;55348.2;53729.8;52153.9;50619.3;49125.1;47670.5;...
    46254.5;44876.2;43534.8;42229.4;40959.1;39723.4;38521.3;37352.2;36215.4;35110.1;...
    34035.7;32991.5;31977;30991.4;30034.3;29104.9;28202.9;27327.6;26478.5;25655.1;...
    24856.9;24083.4;23334.1;22608.7;21906.5;21227.3;20570.6;19936;19323.1;18731.5;...
    18160.9;17610.9;17081.1;16571.2;16080.9;15609.8;15157.7;14724.3;14309.1;13912 ];
elseif strcmp( dataStr, '130607_raw' )
  alphaRadiusVector = [
    333279;326764;323529;260557;196514;153900;124261;96807.6;170481;196651;...
    213216;208154;203402;198932;194722;190749;186998;300449;299276;298071;...
    296839;295579;294295;292985;291655;290305;288937;277143;276363;275565;...
    274743;273900;273040;272163;271266;270357;269439;268505;267557;565984;...
    557358;548983;540856;532973;525315;517869;510636;616664;496753;490087;...
    483602;477286;471121;465123;459274;373547;370707;367916;365175;362478;...
    359829;357223;354655;352129;349646;347201;344794;342417;340082;337781;...
    335511;333277;331070;328898;326754;324642;322558;320497;318468;316465;...
    314489;312536;310607;308704;306825;304967;303134;301320;299527;297758;...
    296009;294280;292566;290876;289204;287549;285914;284294;282693;281107;...
    279540;159994;158753;157531;156328;155142;153972;152818;151682;150561;...
    149455;148365;147289;146228;145181;144149;143129;142122;141130;140150;...
    139182;138227;137283;136351;135431;134522;133625;132738;131861;130996;...
    130141;129296;128459;127634;126818;126010;125213;124424;123644;122873;...
    122110;121356;120609;119871;119142;118420;117705;116997;116298;115606;...
    95210.1;93971.9;92752.6;91552.8;90371.7;89209.2;88064.5;86937;85827.2;84734.4;...
    83657.9;82598;81553.6;80525.2;79511.9;78513.9;77530.7;76561.6;75606.8;74666.2;...
    73739.2;72825.5;71924.9;71037.5;70162.9;69300.4;68450.5;67612.2;66786.2;65971.5;...
    65168.4;64376.5;63595.3;62825.3;62065.9;61316.9;60578.2;59849.4;59130.9;58421.9;...
    57722.7;57032.9;56352.3;55680.9;55018.6;54365.1;53720.4;53084;52456.3;51836.9;...
    51225.5;50622.3;50026.8;49439.4;48859.5;48287.2;47722.4;47164.7;46614.4;46071.3;...
    45535.1;45005.8;44483.2;43967.5;44989.4;44551.5;44116.1;43683.1;43252.7;42824.8;...
    42399.4;41976.4;41555.8;41137.8;40722.2;40309.1;39898.5;39490.3;28531.2;29051.6;...
    29581.3;30120.5;30669.1;31227.5;31795.8;32374.3;32962.9;33561.9;34171.7;34785.7;...
    34409.9;34036.6;33665.9;31323.6;31323.6;31323.3;31322.6;31321.5;31319.9;31317.8;...
    31315.2;31312;31308.3;31303.9;31298.9;31293.2;31286.7;31279.5;31271.5;31262.7;...
    31253.1;31242.6;31231.3;31218.9;31205.7;31191.4;31176.1;31159.9;31142.5;31124.1;...
    31104.5;31083.8;31062;31038.9;31014.7;30989.2;30962.4;30934.3;30904.9;30874.3;...
    30842.2;30808.7;30773.9;30737.7;30700;30660.9;30620.2;30578.1;30534.4;30477.8;...
    28469.3;26572.9;24784.1;23098.8;21512.9;20022.4;18623.9;17313.6;16088.5;14945.1 ];
elseif strcmp( dataStr, '131203_raw' )
  alphaRadiusVector = [
    2.18127e+06;2.49326e+06;2.85033e+06;3.2588e+06;3.72559e+06;4.25848e+06;2.25898e+06;2.47165e+06;2.70924e+06;2.97434e+06;...
    3.26978e+06;1.86267e+06;1.84489e+06;1.827e+06;1.51978e+06;1.66255e+06;1.81931e+06;1.99385e+06;2.18103e+06;2.38172e+06;...
    2.59687e+06;2.8275e+06;3.07475e+06;3.33986e+06;3.62418e+06;3.92919e+06;4.25654e+06;4.60804e+06;4.98565e+06;5.39156e+06;...
    5.8282e+06;6.29824e+06;6.80467e+06;7.35079e+06;4.39642e+06;4.56216e+06;4.71799e+06;4.86345e+06;4.998e+06;5.12147e+06;...
    5.23367e+06;5.33459e+06;7.81825e+06;8.15893e+06;8.1848e+06;7.68951e+06;6.6903e+06;5.89461e+06;5.61156e+06;5.41972e+06;...
    5.23454e+06;5.05581e+06;4.88322e+06;4.71651e+06;4.55543e+06;4.39978e+06;4.24933e+06;4.10387e+06;3.96321e+06;3.82718e+06;...
    1.72492e+06;1.80656e+06;1.3862e+06;1.50142e+06;1.60837e+06;1.7046e+06;1.78795e+06;1.85695e+06;1.91058e+06;1.94848e+06;...
    1.97085e+06;1.97829e+06;1.9719e+06;1.95305e+06;1.92319e+06;1.88401e+06;1.83714e+06;1.7841e+06;1.72638e+06;1.66534e+06;...
    1.6021e+06;1.53771e+06;1.473e+06;1.40871e+06;560979;559725;558572;557534;556586;555757;...
    555023;554400;553874;553458;553140;552935;552827;552832;552935;553153;...
    553464;553893;554423;555074;555822;556696;557663;558761;559961;556825;...
    540398;524813;510016;495945;482556;469796;457630;446014;485033;483145;...
    468419;454466;441238;428676;416742;405386;394576;384271;374444;365059;...
    356092;347516;339309;331447;323913;316685;309749;303085;296682;290522;...
    284595;278887;273389;268087;262974;258039;253275;248671;244222;239919;...
    235756;231727;227825;224044;220381;216828;213382;210038;206792;203639;...
    200577;197600;194706;191892;189154;186489;183895;181369;178909;176511;...
    174174;171895;169673;167504;165389;163324;161309;159339;157417;155537;...
    153701;151906;150151;148434;146755;145112;143505;141931;140390;138881;...
    137024;135430;134457;133848;133518;133402;133446;133605;133847;133882;...
    133459;133049;132651;132264;131888;131524;131170;130827;130493;130171;...
    129858;129555;129261;128976;128701;128434;128176;127926;127684;127451;...
    127225;127007;126797;126594;126399;126210;126029;125855;125688;125527;...
    125373;125225;125084;124948;124819;124696;124579;124467;124360;124260;...
    124166;124076;123992;123914;123840;123772;123709;123650;123458;122886;...
    122312;121736;121161;120582;120006;119429;118853;118278;117704;117130;...
    116558;115987;115419;114852;114288;113725;113165;112607;112053;111498;...
    110949;110403;109860;109319;108782;108247;107717;107190;106667;106146;...
    105630;105117;104607;103015;102404;101810;101231;100666;100116;99577.8;...
    99053.5;98541.1;98041.2;97552.5;97075.5;96608.8;96153.1;95707;95271.2;94843.5;...
    94426.3;94017.8;93618.2;93226.5;92843;92467.3;92099.8;91739.2;91386.2;91039.8;...
    90700.4;90367.3;90041;89720.5;89405.5;89097;88794.5;88497.2;88205.7;87919.1;...
    87638;87361.5;87090.3;86823.4;86561.4;86303.7;86050.7;85801.9;85557;85316.5;...
    85078.9;84846.3;84616.7;84391.5;84169.7;83951.7;83736.9;83525.6;83317.5;83112.7;...
    82910.9;82712.5;82516.7;82324.2;82134.3;81946.9;81762.6;81581.1;81402.2;81226.1;...
    81052.2;80881.2;80712.3;80545.9;80381.8;80220.1;80060.5;79903.3;79748;79595.1;...
    79443.2;79294.4;79147.5;79002.6;78859.5;78718.4;78579.1;78441.7;78306;78172.2;...
    78039.9;77909.6;77780.7;77653.6;77528;77404.1;77281.2;77160.4;77041;76923.2;...
    76806.7;76691.8;76578.1;76466;76355.1;76245.7;76137.4;76030.6;75925;75820.5;...
    75717.3;75615.1;75514.4;75415;75316.6;75219.6;75123.5;75028.6;74934.8;74842.1;...
    74750.4;74659.8;74570.2;74481.8;74394.2;74307.7;74221.8;74137.2;74053.6;73970.9;...
    73889.2;73808.4;73728.4;73649.4;73571.3;73494;73417.5;73341.9;73267.1;73193.2;...
    73120.1;73047.6;72976.1;72905.4;72835.5;72766.4;72698;72630.3;72563.5;72497.3;...
    72431.9;72367.2;72303.1;72239.8;72177.2;72115.3;72054;71993.1;71933.2;71873.9;...
    71815.3;71757.3;71699.9;71643.2;71587.1;71531.4;71476.5;71422.3;71368.5;71315.4;...
    71262.9;71210.9;71159.3;71108.4;71058.2;71008.5;70959.3;70910.7;70862.6;70815 ];
end