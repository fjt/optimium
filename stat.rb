class Array
  def sum
    self.inject(0){|r,e|r+e}
  end

  def mult
    self.inject(1){|r,e|r*e}
  end

  def ave
    self.sum/self.length.to_f
  end

  def variance
    r=self.map{|v|v**2}.ave - self.ave **2
    if r > 0
      r
    else
      0
    end
  end

  def normalise
    sigma=self.std_deviation
    av=self.ave
    self.map{|v|(v-av)/sigma}
  end

  def std_deviation
    Math::sqrt(self.variance)
  end

  def cumdistr
    l=self.length.to_f
    self.sort.reverse.map_with_index{|v, i|[v, (i+1)/l]}
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
    ratio=(step/(self.max-(mi=self.min)))
    self.map{|v|(((mi+v)*ratio).to_i)/ratio-mi}
  end


  def leastsq(isb=t)
    self.depth == 2 and self[0].length == 2 or raise "apply to an array of pairs."
    len=self.length
    mps=self.map{|e|e.mult}.sum
    xs,ys=self.inject([0, 0]){|r,e|r.mplus(e)}
    xss=self.transpose[0].inject(0){|r,e|r+=(e**2)}
    if isb
      [(xs*ys-len*mps)/(xs**2-len*xss).to_f,
	(mps*xs - xss*ys)/(xs**2 - len*xss).to_f] ## returns [slope, shift]
    else
      mps.to_f/xss
    end
  end
end
