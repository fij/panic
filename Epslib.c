// Epslib.c

void EpsInit( FILE *EpsFile, 
	      int boundboxlowleftx, int boundboxlowlefty,
	      int boundboxuprightx, int boundboxuprighty ) {

	fprintf( EpsFile, "%%!PS-Adobe-2.0 EPSF-2.0\n" );
	fprintf( EpsFile, "%%%%BoundingBox: %d %d %d %d\n", 
		 boundboxlowleftx, boundboxlowlefty,
		 boundboxuprightx, boundboxuprighty );
	fprintf( EpsFile, "/cp {closepath} bind def\n" );
	fprintf( EpsFile, "/ef {eofill} bind def\n" );
	fprintf( EpsFile, "/f {fill} bind def\n" );
	fprintf( EpsFile, "/gr {grestore} bind def\n" );
	fprintf( EpsFile, "/gs {gsave} bind def\n" );
	fprintf( EpsFile, "/sa {save} bind def\n" );
	fprintf( EpsFile, "/rs {restore} bind def\n" );
	fprintf( EpsFile, "/l {lineto} bind def\n" );
	fprintf( EpsFile, "/m {moveto} bind def\n" );
	fprintf( EpsFile, "/rm {rmoveto} bind def\n" );
	fprintf( EpsFile, "/n {newpath} bind def\n" );
	fprintf( EpsFile, "/s {stroke} bind def\n" );
	fprintf( EpsFile, "/sh {show} bind def\n" );
	fprintf( EpsFile, "/slc {setlinecap} bind def\n" );
	fprintf( EpsFile, "/slj {setlinejoin} bind def\n" );
	fprintf( EpsFile, "/slw {setlinewidth} bind def\n" );
	fprintf( EpsFile, "/srgb {setrgbcolor} bind def\n" );
	fprintf( EpsFile, "/rot {rotate} bind def\n" );
	fprintf( EpsFile, "/sc {scale} bind def\n" );
	fprintf( EpsFile, "/sd {setdash} bind def\n" );
	fprintf( EpsFile, "/ff {findfont} bind def\n" );
	fprintf( EpsFile, "/sf {setfont} bind def\n" );
	fprintf( EpsFile, "/scf {scalefont} bind def\n" );
	fprintf( EpsFile, "/sw {stringwidth} bind def\n" );
	fprintf( EpsFile, "/tr {translate} bind def\n" );
	fprintf( EpsFile, "/Times-Bold ff 50 scf sf\n" );
}


/* Available eps fonts are: Times-Roman, Times-Italic,
   Times-Bold, Times-Bold-Italic, etc. */

void EpsSetFont( FILE *EpsFile, char *EpsFontStr, double EpsFontHeight ) { 
  fprintf( EpsFile, "/%s ff %g scf sf\n", EpsFontStr, EpsFontHeight ); }



void EpsDrawString( FILE *EpsFile, double deg, double x, double y, char *TextStr ) {
  /* draws rotated text onto picture -- used until now only as last
     command in file 
     deg: angel of rotation in degrees */

  fprintf( EpsFile, "%.2f rot\n", deg );
  fprintf( EpsFile, "n %g %g m (%s) sh s\n", 
	   /* formula works for deg=0.0, 
	      for other angles, try -- formula may not work */
	   x*cos(deg)+y*sin(deg),
	   y*cos(deg)-x*sin(deg), 
	   TextStr ); 
}


void EpsDrawRectangle( FILE *EpsFile, 
		       double lowleftx, double lowlefty,
		       double uprightx, double uprighty ) {

  fprintf( EpsFile, "n %g %g m\n", lowleftx, lowlefty );
  fprintf( EpsFile, "%g %g l\n", lowleftx, uprighty );
  fprintf( EpsFile, "%g %g l\n", uprightx, uprighty );
  fprintf( EpsFile, "%g %g l\n", uprightx, lowlefty );
  fprintf( EpsFile, "cp s\n" );
}


void EpsFillRectangle( FILE *EpsFile, 
		       double lowleftx, double lowlefty,
		       double uprightx, double uprighty ) {

  fprintf( EpsFile, "n %g %g m\n", lowleftx, lowlefty );
  fprintf( EpsFile, "%g %g l\n", lowleftx, uprighty );
  fprintf( EpsFile, "%g %g l\n", uprightx, uprighty );
  fprintf( EpsFile, "%g %g l\n", uprightx, lowlefty );
  fprintf( EpsFile, "cp f s\n" );
}



void EpsDrawCircle( FILE *EpsFile, double centerx, double centery, double rad )
{
  fprintf( EpsFile, "n %g %g m %g %g %g 0 360 arc s\n", \
	   centerx + rad, centery, centerx, centery, rad );
}


void EpsFillCircle( FILE *EpsFile, double centerx, double centery,
		    double rad ) {
  fprintf( EpsFile, "n %g %g m %g %g %g 0 360 arc f s\n", \
	   centerx + rad, centery, centerx, centery, rad );
}


void EpsSetLinewidth( FILE *EpsFile, double linewidth ) {
  fprintf( EpsFile, "%g slw\n", linewidth ); }




void EpsSetRgb( FILE *EpsFile, double red, double green, double blue ) {
  fprintf( EpsFile, "%g %g %g srgb\n", red/65535.0, green/65535.0, blue/65535.0 ); }

