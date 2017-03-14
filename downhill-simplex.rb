class Array
  def mplus (arg)
    if self[0].is_a? Numeric
      self.map.with_index{|v, i| v + arg[i]}
    else
      self.map.with_index{|v, i|v.mplus(arg[i])}
    end
  end

  def abs
    if self.first.class == self.class
      sqrt(self.map{|e|
             e.abs**2}.sum)
    else
      sqrt(self.map{|v|v**2}.sum)
    end
  end

  def sum
    self.inject(0){|r,v|v+=r}
  end

  def ave
    self.sum.to_f/self.length
  end

  def dot(arg)
    if arg.is_a? Numeric
      if self.first.is_a? Numeric
        self.map{|v| v * arg}
      else
        self.map{|e|e.dot(arg)}
      end
    else
      if arg.first.is_a? Numeric ## arg is vector
	if self.first.is_a? Numeric ## self is vector
	  self.inject_with_index(0){|ret, v, i|ret += v * arg[i]}
	else ##  matrix dot vector
	  self.map{|row|row.dot(arg)}
#	  (NMatrix.to_na(self) * NVector.to_na(arg)).to_a
	end
      else ## matrix dot matrix
	self.map{|v|arg.transpose.map{|u| v * u}}
#	(NMatrix.to_na(self) * NMatrix.to_na(arg)).to_a
      end
    end
  end
end

class Proc
  def dhsmplx(simplex, count=1000, hist=[], lt={})
    hist.push(simplex)
    ss=simplex[0]; size= simplex.map{|v|ss.mplus(v.dot(-1)).abs}.inject(0){|r, v|r+=v}
    count-=1
    if count == 0
      hist
    else
      ref = simplex.map{|v|
	[ (lt[v] or lt[v] = self.call(*v)), v]
      }.sort
      best = ref[0]; worse = ref[-2]; worst = ref[-1]
      nsimplex = ref[0..-2].map{|e|e[1]} ## remove worst
      cp = nsimplex.transpose.map{|ax|ax.ave} ## centroid
      cr = (cp.dot(2.0)).mplus(worst[1].dot(-1)) ## reflection of the worst
      rv = (lt[cr] = self.call(*cr))
      if best[0] <= rv and rv <= worse[0] ## case intermediate
        dhsmplx(nsimplex.push(cr), count, hist, lt)
      else
        if rv < best[0] ## best
          ep = (cr.dot(2.0)).mplus(cp.dot(-1)) ## expand
          if (lt[ep] = self.call(*ep)) < rv
            nsimplex.unshift(ep)
          else
            nsimplex.unshift(cr)
          end
          dhsmplx(nsimplex, count, hist, lt)
        else ## case worst
            ctrp = (cp.mplus(worst[1])).dot(0.5) ## contract point
          if (lt[ctrp] = self.call(*ctrp)) <= worst[0]
            dhsmplx(nsimplex.push(ctrp), count, hist, lt)
          else
            nsimplex = [best[1]]+(ref[1..-1].map{|e|e[1]}.map{|v|(best[1].mplus(v)).dot(0.5)})
            dhsmplx(nsimplex, count, hist, lt)
          end
        end
      end
    end
  end
end
