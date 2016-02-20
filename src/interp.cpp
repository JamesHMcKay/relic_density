#include "interp.hpp"


int Base_interp::locate(const double x)
{
  int ju, jm, jl;
  if (n<2 || mm<2 || mm>n) throw("locate size error");

  bool ascnd = (xx[n-1] >= xx[0]);
  jl=0;
  ju=n-1;

  while (ju-jl >1 )
  {
  jm= (ju+jl) >>1;
//  cout << "middle point = " << jm << endl;
//  cout << "upper limit = " << ju << endl;
//  cout << "lower limit = " << jl << endl;
      if ((x >= xx[jm]) == ascnd){
        jl=jm;
       }
      else{
        ju=jm;
        }
  }

  cor = abs (jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return MAX(0,MIN(n-mm, jl-((mm-2)>>1)));
}

int Base_interp::hunt(const double x)
{
int jl=jsav, jm, ju, inc=1;
if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
bool ascnd=(xx[n-1] >= xx[0]);
if (jl < 0 || jl > n-1)
{
jl=0;
ju=n-1;
}
else
{
  if ((x >= xx[jl]) == ascnd)
  {
    for (;;)
    {
     ju = jl + inc;
     if (ju >= n-1) { ju = n-1; break;}
     else if ((x < xx[ju]) == ascnd) break;

     else
     {
      jl = ju;
      inc += inc;
     }
    }
   }
   else
   {
    ju = jl;
    for (;;)
    {
      jl = jl - inc;
      if (jl <= 0) { jl = 0; break;}
      else if ((x >= xx[jl]) == ascnd) break;
      else
      {
        ju = jl;
        inc += inc;
      }
    }
   }
}
while (ju-jl > 1)
{
    jm = (ju+jl) >> 1;
    if ((x >= xx[jm]) == ascnd)
    jl=jm; else
    ju=jm;
}
cor = abs(jl-jsav) > dj ? 0 : 1;
jsav = jl;
return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}


double Poly_interp::rawinterp(int jl, double x)
{
int i,m,ns=0;
double y,den,dif,dift,ho,hp,w;
const double *xa = &xx[jl], *ya = &yy[jl]; std::vector<double> c(mm),d(mm);
dif=abs(x-xa[0]);
for (i=0;i<mm;i++)
{
if ((dift=abs(x-xa[i])) < dif)
{
ns=i;
dif=dift;
}
c[i]=ya[i];
d[i]=ya[i];
}
y=ya[ns--];
for (m=1;m<mm;m++)
{

for (i=0;i<mm-m;i++)
{
    ho=xa[i]-x;
    hp=xa[i+m]-x;
    w=c[i+1]-d[i];
    if ((den=ho-hp) == 0.0) throw("Poly_interp error");
    den=w/den;
    d[i]=hp*den;
    c[i]=ho*den;
}
y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
}
return y;
}







