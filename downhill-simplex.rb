class Array
  def mplus (arg)
    if self[0].is_a? Numeric
      self.map_with_index{|v, i| v + arg[i]}
    else
      self.map_with_index{|v, i|v.mplus(arg[i])}
    end
  end    

  def *(arg)
    if arg.is_a? Numeric
      self.map{|v| v * arg}
    else
      if arg[0].is_a? Numeric
	if self[0].is_a? Numeric
	  self.inject_with_index(0){|ret, v, i|ret += v * arg[i]}
	else
	  self.map{|row|row * arg}
#	  (NMatrix.to_na(self) * NVector.to_na(arg)).to_a
	end
      else
	self.map{|v|arg.transpose.map{|u| v * u}}
#	(NMatrix.to_na(self) * NMatrix.to_na(arg)).to_a
      end
    end
  end
end



class Proc
  def dhsmplx(simplex, hist=[], lt={})
    hist.push(simplex)
    ss=simplex[0]; size= simplex.map{|v|ss.mplus(v*-1).abs}.inject(0){|r, v|r+=v}
    if size < 0.00000001
      hist
    else
      #      ref = simplex.map{|v|[self.call(*v), v]}.sort
      ref = simplex.map{|v|
	[ (lt[v] or lt[v] = self.call(*v)), v]
      }.sort
      best = ref[0]; worse = ref[-2]; worst = ref[-1]
      nsimplex = ref[0..-2].map{|e|e[1]} ## remove worst
      cp = simplex.transpose.map{|ax|ax.ave} ## centroid
      cr = (cp*2.0).mplus(worst[1]*(-1)) ## reflection of the worst
      rv = self.call(*cr)
      if best[0] <= rv and rv < worse[0] ## case intermediate
        dhsmplx(nsimplex.push(cr), hist, lt)
      else
        if rv < best[0] ## best
          ep = (cp*3.0).mplus(worst[1]*(-2))
          if self.call(*ep) < best[0]
            nsimplex.push(ep)
          else
            nsimplex.push(cr)
          end
          dhsmplx(nsimplex, hist, lt)
        else ## case worst
            ctrp = (cp.mplus(worst[1]))*0.5 ## contract point
          if self.call(*ctrp) < worse[0]
            dhsmplx(nsimplex.push(ctrp), hist, lt)
          else
            nsimplex = ref[1..-1].map{|e|e[1]}.map{|v|(best[1].mplus(v))*0.5}.push(best[1]) ## conctact others to the best point
            dhsmplx(nsimplex, hist, lt)
          end
        end
      end
    end
  end
end
