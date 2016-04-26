module RandomVariate

  class RandomVariateGenerator
    def new(sample)
      deriv=sample.each_cons(2).map{|e|e.last - e.first}
      dderiv=deriv.each_cons(2).map{|e|e.last - e.first}
    end
  end

  def gamma(opt={})
    scale = (opt[:scale] or 1.0)
    shape = (opt[:shape] or 2.0)

    sg_magicconst = 1.0 + Math::log(4.5)

    raise "Gamma params should be positive, set to 1.0, 2.0 by default." if scale <= 0 or shape <= 0
    if shape > 1.0
      ainv = Math::sqrt(shape*2.0 - 1.0)
      b = shape - Math::log(4)
      c = shape + ainv

      while true
	if 1.0e-7 < (u1 = rand(0)) and u1 < 0.9999999
	  u2 = 1.0 - rand(0)
	  v = Math::log(u1/(1.0-u1))/ainv
	  x = shape*Math::E**(v)
	  z = u1*u1*u2
	  r = b + c*v - x
	  return(x * scale) if r + sg_magicconst - 4.5*z >= 0.0 or r >= Math::log(z)
	end
      end

    else if shape == 1.0
	   u = rand(0)
	   while u = rand(0) <= 1.0e-7
	   end
	   return (-Math::log(u) * scale)
	 else
	   while true
	     u = rand(0)
	     b = (Math::E + shape)/Math::E
	     p = b*u
	     if p <= 1.0
	       x = p**(1.0/shape)
	     else
	       x = Math::log((b-p)/shape)
	     end
	     u1 = rand(0)
	     if p > 1.0
	       if u1 <= x**(shape -1.0)
		 break
	       end
	     else if u1 <= E**(-x)
		    break
		  end
	     end
	     return ( x * beta)
	   end
	 end
    end
  end


  def cdfinvert(&cdfinv)
    cdfinv.call(rand(0))
  end

  def power(opt={})
    ## gives variates obeying given power law distribution.
    ## the lower bound and the power can be given by option, which are set to 1.0 and 2.0 by default.
    nom = (opt[:nom] or 1.0)
    power = (opt[:power] or 2.0)
    cdfinvert{|x|nom/(1-x)**(1.0/(power-1))}
  end

  def clgauss(mean=0, var=1.0)
    sigma=Math::sqrt(var)
    sigma*(Array.new(12){rand(0)}.inject(0){|r, v|r+=v}-6.0) + mean
  end

  def bmgauss(arg=nil)
    alpha, beta = rand(0), rand(0)
    if arg.nil?
      Math::sqrt(Math::log(alpha)*(-2))*Math::cos(2*Math::PI*beta)
    else
      var = (arg[:var] or 1)
      mean = (arg[:mean] or 0)
      if arg[:dim] == 2
	[Math::sqrt(Math::log(alpha)*(-2))*Math::cos(2*Math::PI*beta), Math::sqrt(Math::log(alpha)*(-2))*Math::sin(2*Math::PI*beta)].map{|v|v*sigma + mu}
      else
	Math::sqrt(Math::log(alpha)*(-2))*Math::cos(2*Math::PI*beta) * var + mean
      end
    end
  end

  alias :gauss :bmgauss ## use Box-Muller method
  alias :gaus :bmgauss


  def erf(x)
    ## based on http://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/#
    a1, a2, a3, a4, a5, p =  0.254829592,\
    -0.284496736,\
    1.421413741,\
    -1.453152027,\
    1.061405429,\
    0.3275911

    sign = x > 0
    x = x.abs

    tee = 1.0/(1.0 + p*x)
    wai = 1.0 - ((((a5*tee + a4)*tee + a3)*tee + a2)*tee + a1)*tee*(E**(-x*x))
    if sign; wai; else; -wai;end
  end

  def erfinv(p)
    raise  "arg must be positive and less than 1.0" if p >= 1 or p <= 0
    a_1 = -3.969683028665376e+01
    a_2 =  2.209460984245205e+02
    a_3 = -2.759285104469687e+02
    a_4 =  1.383577518672690e+02
    a_5 = -3.066479806614716e+01
    a_6 =  2.506628277459239e+00

    b_1 = -5.447609879822406e+01
    b_2 =  1.615858368580409e+02
    b_3 = -1.556989798598866e+02
    b_4 =  6.680131188771972e+01
    b_5 = -1.328068155288572e+01

    c_1 = -7.784894002430293e-03
    c_2 = -3.223964580411365e-01
    c_3 = -2.400758277161838e+00
    c_4 = -2.549732539343734e+00
    c_5 =  4.374664141464968e+00
    c_6 =  2.938163982698783e+00

    d_1 =  7.784695709041462e-03
    d_2 =  3.224671290700398e-01
    d_3 =  2.445134137142996e+00
    d_4 =  3.754408661907416e+00

    p_low  = 0.02425
    p_high = 1 - p_low

    case p
    when (0..p_low)
      q = sqrt(-2*log(p))
      x = (((((c_1*q+c_2)*q+c_3)*q+c_4)*q+c_5)*q+c_6) /\
      ((((d_1*q+d_2)*q+d_3)*q+d_4)*q+1)
    when (p_low..p_high)
      q = p - 0.5
      r = q*q
      x = (((((a_1*r+a_2)*r+a_3)*r+a_4)*r+a_5)*r+a_6)*q /
           (((((b_1*r+b_2)*r+b_3)*r+b_4)*r+b_5)*r+1)

  when (p_high..1.0)
      q = sqrt(-2*log(1-p))
      x = -(((((c_1*q+c_2)*q+c_3)*q+c_4)*q+c_5)*q+c_6) /
      ((((d_1*q+d_2)*q+d_3)*q+d_4)*q+1)
    end
  end
end ## end module Random

class Numeric
  def ub(b)
    if self > b
      b
    else
      self
    end
  end

  def lb(b)
    if self < b
      b
    else
      self
    end
  end
end

class Array

  def ent(base=2)
    len=self.length.to_f
    - self.inject({}){|h, e|if h[e]; h[e]+=1; else; h[e]=1; end; h}.to_a.transpose.last.map{|v|(p=v/len)*Math::log(p,base)}.sum
  end

  def distdist(arg)
    if arg.length > self.length
      arg=arg.sample(self.length).sort
      base=self.sort
    else
      arg=arg.sort
      base=self.sample(arg.length).sort
    end
    len=base.length.to_f
    base.each_with_index.inject(0){|r, e|
      i=e.last
      val=e.first
      r+(i - arg.filter{|v|v< val}.length).abs/len}/len
  end

  def ksdist(arg)
    if arg.length > self.length
      arg=arg.sample(self.length).sort
      base=self.sort
    else
      arg=arg.sort
      base=self.sample(arg.length).sort
    end
    base.map.with_index{|val, i|
      (i - arg.filter{|v|v< val}.length).abs}.max/base.length.to_f
  end

  def rmap(&bl)
    self.map{|e|
      if e.class == self.class
	e.rmap(&bl)
      else
        bl.call(e)
      end}
  end

  def filter(&bl)
     self.inject([]){|r, e| if e=bl.call(e)
                              r.push(e)
                            else
                              r
                            end}
#    self.map(&bl).map{|e|nil if e == false}.compact
  end


  def abs
    Math::sqrt(self.inject(0){|r, v|r+= v*v})
  end

  def eig_pow
    eps=0.0001
    l=0.0; ol=1.0
    if self[0].class == self.class
##      x=[1];(self[0].length - 1).times{x.push(0)}
      x=self.map{|r|r.sum}
      while (l - ol).abs/ol > eps
        ol=l
        y= self.prod(x)
        l=y.abs / x.abs
        x=y/l
      end
      xa=x.abs
      [x.map{|v|v/xa}, l]
    else
      x=self.map.with_index{|r, i|r.length}.to_cv
      while (l - ol).abs/ol > eps
        ol=l
        y= self.cmprod(x)
        xabs, yabs = [x, y].map{|v|v.to_a.transpose[1].abs}
        l=yabs / xabs
        x= y.cmap{|k, v| v/l}
      end
      xabs=0; x.each{|k, v|xabs+=v}
      [x.cmap{|k, v|v/xabs}, l]
    end
  end

  def prod(arg, &block)
    if arg[0].is_a? Numeric
      if self[0].is_a? Numeric
#        self.inject_with_index(0){|ret, v, i|ret+= block.call(v, arg[i])}
        [self, arg].transpose.inject(0){|r, e|r+e.first*e.last}
      else
        self.map{|r| r.prod(arg, &block)}
      end
    else
      self.map{|v|arg.map{|vv| v.prod(vv, &block)}}
    end
  end



  def sum
    self.inject(0){|r,e|r+e}
  end

  def mult
    self.inject(1){|r,e|r*e}
  end

  def ave
    self.sum/self.length.to_f
  end

  def movave(len)
    len=len.to_f
    sm=self[0..(len-1)].sum
    ret=[sm/len]
    (0..(self.length - len -1)).each{|i|
      sm=sm-self[i]+self[i+len]
      ret.push(sm/len)}
    ret
  end

  def rpick
    self[rand(self.length)]
  end

  def tocopula
    en=self.length.to_f
    self.transpose.map{|row|
      i=0
      im = row.sort.inject({}){|r,e|i+=1; r[e]=i/en; r}
      row.map{|v|im[v]}}.transpose
  end

  def leastsq(isb=true)
    self.depth == 2 and self[0].length == 2 or raise "apply to an array of pairs."
    len=self.length
    mps=self.map{|e|e.mult}.sum
    xs,ys=self.transpose.map{|e|e.sum}
    xss=self.transpose[0].inject(0){|r,e|r+=(e**2)}
    if isb
      [(xs*ys-len*mps)/(xs**2-len*xss).to_f,
	(mps*xs - xss*ys)/(xs**2 - len*xss).to_f] ## returns [slope, shift]
    else
      mps.to_f/xss
    end
  end

#   def variance_obs
#     r=self.map{|v|v**2}.ave - self.ave **2
#     if r > 0
#       r
#     else
#       0
#     end
#   end

  def variance
    dnm = (self.length - 1).to_f
    av=self.ave
    (self.inject(0){|r, v|r + (v - av)**2})/ dnm
  end

  def std_deviation
    Math::sqrt(self.variance)
  end

  def normalise
    sigma=self.std_deviation
    av=self.ave
    self.map{|v|(v-av)/sigma}
  end

  alias :normalize :normalise

  def max(&block)
    block=lambda{|a, b| a > b} if block.nil?
    self.reduce{|i, j|j = i if j.nil?
      if block.call(i, j)
	i
      else
	j
      end}
  end

  def min
    self.reduce{|i, j|j = i if j.nil?
      if i < j
	i
      else
	j
      end}
  end

  def cumdistr(dir=:increase)
    l=self.length.to_f
    if dir == :increase
      self.sort.map_with_index{|v, i|[v, (i+1)/l]}
    else
      self.sort.reverse.map_with_index{|v, i|[v, (i+1)/l]}
    end
  end

  def distr
    (l=self.sort).inject([[l[0], 0]]){|r, v|
      if v == r[-1][0]
	r[-1][1]+=1
      else
	r.push([v, 1])
      end
      r}
  end

  def fold(step=100)
    ratio=(step.to_f/(self.max-(mi=self.min)))
    self.map{|v|(((mi+v)*ratio).to_i)/ratio-mi}
  end

  def mle_power(xmin = self.min, error = false)
    ## calculate cumulative distribution power, not density like my random variate function.
    ## Hill's estimate
    xmin= xmin.to_f
    d= self.map{|e|e if e >= xmin}.compact
    l= d.length
    pow = l / d.inject(0){|r, i|r+= Math::log(i/xmin)}
    if not error
      pow
    else
      error = Math::sqrt(l + 1) / d.inject(0){|r, i|r+=Math::log(i/xmin)}
      [pow, [pow - error .. pow + error]]
    end
  end


  def pearson_corl
    self.transpose.map{|row|row.normalize}.transpose.inject(0){|p, row|p+row.mult} / \
    (self.length - 1).to_f
  end

  def value2rank(sorted=nil)
    self.map{|v|
      if v.is_a? Enumerable
        v.value2rank(sorted)
      else
        sorted = self.sort unless sorted
        sorted.index(v)
      end}
  end

  def v2rank
    p=nil; pl=0; ret=Hash.new
    self.sort.map.with_index{|v, i|
      if p
        if p < v
          p=v; pl=i; ret[v]=i
        else
          ret[v]=pl
        end
      else
        p=v; ret[v]=0
      end
    }
    ret
  end

  def spearman_corl
     ecs= self.map{|e|e[0]}.v2rank
     wai= self.map{|e|e[1]}.v2rank
     self.map{|v|[ecs[v[0]], wai[v[1]]]}.pearson_corl
  end

  def revpairs_length
    if self.length > 1
      pivot=self[0]
      head, tail = [], []
      t_h_revs=0
      self[1..-1].each{|e|
      if e > pivot
	  tail.push(e)
	else
	  head.push(e)
	  t_h_revs+=tail.length
	end
      }
      t_h_revs + head.inject(0){|r, e|if e == pivot; r+0.5;else;r+1;end} + head.revpairs_length + tail.revpairs_length
    else
      0
    end
  end

  def concordant_pairs
    if self.length > 1
      pivot=self[0]
      head, tail = [], []
      t_h_cons=0
      self[1..-1].each{|e|
      if e > pivot
        tail.push(e)
        t_h_cons+=head.length
      else 
        head.push(e)
      end}
      t_h_cons + tail.length + head.concordant_pairs + tail.concordant_pairs
    else
      0
    end
  end

  def discordant_pairs
    if self.length > 1
      pivot=self[0]
      head, tail = [], []
      t_h_revs=0
      self[1..-1].each{|e|
      if e < pivot
        head.push(e)
        t_h_revs+=tail.length
      else
        tail.push(e)
      end
      }
      t_h_revs + head.length + head.discordant_pairs + tail.discordant_pairs
    else
      0
    end
  end

  def tie_pairs
    if self.length > 1
      piv=self.first
      tp=0
      hd, tl = [], []
      self[1..-1].each{|e|
        if e > piv
          tl.push(e)
        else if e < piv
               hd.push(e)
             else
               tp+=1
             end
        end
      }
      tp + tl.tie_pairs + hd.tie_pairs
    else
      0
    end
  end

  def kendall_tau
    nomin = self.transpose.map{|e|el=e.length; (el-1)*el/2 - e.tie_pairs}.inject(1){|r,e|r*e}
    ary=self.sort.transpose.last
    (ary.concordant_pairs - ary.discordant_pairs) / sqrt(nomin)
  end

  def kendall_corl
    pairs= (en=self.length)*(en-1)/2.0
    (pairs - 2*self.sort.transpose[1].revpairs_length)/pairs
  end

  def jaccard(other)
    (self & other).uniq.length.to_f/(self + other).uniq.length
  end

  def mode
    self.inject({self.first => 1}){|h,e|if h[e]; h[e]+=1; else; h[e]=1;end;h}.to_a.sort{|a,b|a.last <=> b.last}.last
  end

  def entropy_continuous(base=2.0)
    d=self.sort.each_cons(2).map{|e|1/(e.last - e.first) if e.last> e.first}.compact
    s=d.sum.to_f
    d.inject(0){|h, v|h+=(v/s)*log(s/v, base)}
  end

end
