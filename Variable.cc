#include "Variable.hh"
#include <cmath> 

Variable::Variable (std::string n) 
  : Indexable(n) 
  , numbins(100)
  , fixed(false)
  , blind(0)
  , error(0.0)
  , error_pos(0.0)
  , lowerlimit(0.0)
  , upperlimit(0.0)
  , error_neg(0.0)
  , gcc(0.0)
{
} 

Variable::Variable (std::string n, fptype v) 
  : Indexable(n, v)
  , error(0.002) 
  , error_pos(0.0)
  , error_neg(0.0)
  , gcc(0.0)
  , lowerlimit(v - 0.01)
  , upperlimit(v + 0.01)
  , numbins(100)
  , fixed(true)
  , blind(0)
{
}

Variable::Variable (std::string n, fptype dn, fptype up) 
  : Indexable(n)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
  , error_pos(0.0)
  , error_neg(0.0)
  , gcc(0.0)
  , error(0.0)
{
}

Variable::Variable (std::string n, fptype v, fptype dn, fptype up) 
  : Indexable(n, v)
  , error(0.1*(up-dn))
  , error_pos(0.0)
  , error_neg(0.0)
  , gcc(0.0)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
{
}

Variable::Variable (std::string n, fptype v, fptype e, fptype dn, fptype up) 
  : Indexable(n, v)
  , error(e)
  , error_pos(0.0)
  , error_neg(0.0)
  , gcc(0.0)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
{
}

Variable::~Variable () {
}

std::ostream& operator<<(std::ostream &oss, const Variable &var)
{
  oss << "Variable with:" << std::endl
  << "\tname = " << var.name << std::endl
  << "\tvalue = " << var.value << std::endl
  << "\terror = " << var.error << std::endl
  << "\terror_neg = " << var.error_neg << std::endl
  << "\terror_pos = " << var.error_pos << std::endl
  << "\tlowerlimit = " << var.lowerlimit << std::endl
  << "\tupperlimit = " << var.upperlimit << std::endl
  << "\tnumbins = " << var.numbins << std::endl
  << "\tindex = " << var.index << std::endl
  << "\tfixed = " << (var.fixed ? "true" : "false") << std::endl;
  return oss;
}

