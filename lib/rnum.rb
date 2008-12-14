#
#  rnum.rb
#  Ruby Numerical Library - RNum
#    (C) Copyright 2006- by Olaf Trygve BERGLIHN
#
#  This program is free software.  You can distribute/modify this program
#  under the same terms as Ruby itself.  There is absolutely no warranty with
#  this software.
#
# = rnum.rb

#
# The RNum class is the base class for Ruby Numeric Library.  It implements
# constructors and methods shared by vector and matrix subclasses.
#
# TODO: Ruby classes for complex case needs to be written.
#
module RNum
    require 'ratlas'
    
    #
    # Creates a real matrix given a array of rows:
    #   RNum::matrix([[1, 2, 3], [4, 5, 6]])
    #   => RNum::Matrix[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    #
    def matrix(arg); return Matrix::new(arg); end

    # Creates a real vector given an array:
    #   RNum::vector([1, 2, 3, 4])
    #   => RNum::Vector[1.0, 2.0, 3.0, 4.0]
    # 
    def vector(arg); return Vector::new(arg); end
    
    #
    # Creates a matrix with all elements set to 1.0.  Argument is either a
    # Fixnum or an array of 2 Fixnums.
    #   RNum::ones(2,2)
    #   => RNum::Matrix[[1.0, 1.0], [1.0, 1.0]]
    #   RNum::ones(3,2)
    #   => RNum::Matrix[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
    #   
    def ones(m, n=nil)
        arg = n.nil? ? m : [m,n]
        case arg
        when Array
            Matrix::ones(arg)
        when Fixnum
            Vector::ones(arg)
        else
            raise ArgumentError
        end
    end
     
    # 
    # Creates a matrix with all elements set to 0.0.  Arguments as in
    # ones.
    # 
    def zeros(m, n=nil)
        arg = n.nil? ? m : [m,n]
        case arg
        when Array
            Matrix::zeros(arg)
        when Fixnum
            Vector::zeros(arg)
        else
            raise ArgumentError
        end
    end
    
    # Creates a matrix with 1.0 on the diagonal.  Argumens as in ones.
    # 
    def eye(m, n=nil)
        arg = n.nil? ? m : [m,n]
        case arg
        when Array
                   Matrix::eye(arg)
               when Fixnum
               Vector::eye(arg)
               end
    end
#          def RNum::cmatrix(arg); return self::Cmatrix[arg]; end
#          def RNum::cvector(arg); return self::Cvector[arg]; end
#          def RNum::cones(arg); return self::Cmatrix.con
    
    
    #
    # Allocates memory for a matrix or vector.  Returns a real matrix if
    # size is a two elemene array of Fixnum, returns a real vector if size
    # is a Fixnum.  Does not set memory to zero.
    #
    def alloc(m, n=nil)
        arg = n.nil? ? m : [m,n]
        case arg
            when Fixnum
                stor = RAtlas::alloc([arg,1])
                return Vector.alloc(stor)
            when Array
                stor = RAtlas::alloc(arg)
                return Matrix.alloc(stor)
            else
                raise ArgumentError
        end
    end

    #
    # Concatenates vectors [self, *args]
    #
    def veccat(*args)
        if args.size < 2
          case args.first
          when Numeric
            return vector(args)
          when Vector
            return vector(args[0])
          else
            raise ArgumentError, "Was given only one argument and expected Numeric or Vector"
          end
        end
        dim = args.inject(0){|d, a|
            case a
            when Numeric
                d += 1
            when Vector
                d += a.length
            else
                raise ArgumentError, "Expect Numeric or Vector."
            end
        }
        v = alloc(dim)
        args.inject(0){|i, a|
            case a
            when Numeric
                v[i] = a
                i += 1
            when Vector
                v[i..i+a.length-1] = a
                i += a.length
            end
        }
        return v
    end
    
    #
    # Concatenates horizontally, i.e. stacks elements horizontally.
    #
    def hcat(*args)
        args = args.flatten
        dims = args.inject([0,0]){|d, a|
            case a
            when Numeric
                if d[0] > 1
                    raise ArgumentError, 
                        "Conflicting number of rows."
                end
                d[0] = 1
                d[1] += 1
            when Vector
                if d[0] > a.length
                    raise ArgumentError,
                        "Conflicting number of rows."
                end
                d[0] = a.length    
                d[1] += 1
            when Matrix
                if d[0] > a.size[0]
                    raise ArgumentError,
                        "Conflicting number of rows."
                end
                d[0] = a.size[0]
                d[1] += a.size[1]
            end
            d
        }
        m = alloc(dims)
        colidx = 0
        args.inject(0){|c, a|
            case a
            when Numeric
                if dims[0] != 1
                    raise ArgumentError,
                        "Conflicting number of rows."
                end
                m[0, c] = a
                c += 1
            when Vector
                if dims[0] != a.length
                    raise ArgumentError,
                        "Conflicting number of rows."
                end
                m[0..-1,c] = a
                c += 1
            when Matrix
                if dims[0] != a.size[0]
                    raise ArgumentError,
                        "Conflicting number of rows."
                end
                m[0..-1,c..c+a.size[1]-1] = a
                c += a.size[1]
            end
        }
        return m
    end

    #
    # Concatenates vertically, i.e. stacks arguments vertically.
    #
    def vcat(*args)
        return hcat(*args) if args.size < 2
        dims = args.inject([0,0]){|d, a|
            case a
            when Numeric
                if d[1] > 1
                    raise ArgumentError, 
                        "Conflicting number of columns."
                end
                d[1] = 1
                d[0] += 1
            when Vector
                if d[1] > 1
                    raise ArgumentError,
                        "Conflicting number of columns."
                end
                d[1] = 1    
                d[0] += a.length
            when Matrix
                if d[1] > a.size[1]
                    raise ArgumentError,
                        "Conflicting number of columns."
                end
                d[1] = a.size[1]
                d[0] += a.size[0]
            end
            d
        }
        m = alloc(dims)
        args.inject(0){|r, a|
            case a
            when Numeric
                if dims[1] != 1
                    raise ArgumentError,
                        "Conflicting number of columns."
                end
                m[r, 0] = a
                r += 1
            when Vector
                if dims[1] != 1
                    raise ArgumentError,
                        "Conflicting number of columns."
                end
                m[r..r+a.length-1, 0] = a
                r += a.length
            when Matrix
                if dims[1] != a.size[1]
                    raise ArgumentError,
                        "Conflicting number of columns."
                end
                m[r..r+a.size[0]-1,0..-1] = a
                r += a.size[0]
            end
        }
        return m
    end

    #
    # Concatenates elements.
    #   cat([d,d],a]) = vcat(hcat(d,d),a)
    #
    def cat(*args); vcat(*args.collect{|r| hcat(r)}); end

    #
    # Helper class method for unary operators.
    #
    def RNum::unop(arg, meth)
        case arg
        when RNum::Base
            arg.map{|i| meth.call(i)}
        when Numeric
            meth.call(arg)
        else
            raise ArgumentError
        end
    end
    def sin(arg);   RNum::unop(arg, Math.method(:sin)); end
    def cos(arg);   RNum::unop(arg, Math.method(:cos)); end
    def tan(arg);   RNum::unop(arg, Math.method(:tan)); end
    def asin(arg);  RNum::unop(arg, Math.method(:asin)); end
    def acos(arg);  RNum::unop(arg, Math.method(:acos)); end
    def atan(arg);  RNum::unop(arg, Math.method(:atan)); end
    def sinh(arg);  RNum::unop(arg, Math.method(:sinh)); end
    def cosh(arg);  RNum::unop(arg, Math.method(:cosh)); end
    def tanh(arg);  RNum::unop(arg, Math.method(:tanh)); end
    def asinh(arg); RNum::unop(arg, Math.method(:asinh)); end
    def acosh(arg); RNum::unop(arg, Math.method(:acosh)); end
    def atanh(arg); RNum::unop(arg, Math.method(:atanh)); end
    def log10(arg); RNum::unop(arg, Math.method(:log10)); end
    def exp(arg);   RNum::unop(arg, Math.method(:exp)); end
    def sqrt(arg)
           case arg
           when RNum::Base
               arg.clone.sqrt!
           when Numeric
               Math::sqrt(arg)
           else
               raise ArgumentError
           end
    end
    def log(arg)
        case arg
        when RNum::Base
            arg.clone.log!
        when Numeric
            Math::log(arg)
        else
            raise ArgumentError
        end
    end
    def compare(left, op, right)
      unless left.class == right.class
        raise ArgumentError, "Can not compare unequal RNum types."
      end
      case left
      when RNum::Vector
        left.zip(right){|_a,_b| _a.send(op,_b) ? 1 : 0 }.to_a.map!{|e|
          e > 0.0 ? true : false 
        }
      else 
        left.zip(right){|_a,_b| _a.send(op,_b) ? 1 : 0 }.to_a.map!{|e|
          e.map!{|ee| e > 0.0 ? true : false }
        }
      end
    end
    module_function :matrix, :vector, :ones, :zeros, :eye, :alloc, :sin,
               :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh,
        :acosh, :atanh, :log10, :exp, :sqrt, :log, :veccat, :hcat, 
        :vcat, :cat, :compare

    class Base
        attr_reader :storage
        attr_writer :storage
        private_class_method :new

        # 
        # Return memory address of data storage.
        #
        def memadr; RAtlas::memadr(@storage);end

        #
        # Return memory size of data storage.
        def memsize; RAtlas::memsize(@storage);end

        #
        # Helper method for unary operators
        #
        def unop!(meth); map!{|i| meth.call(i)}; end
        def unop(meth); map{|i| meth.call(i)}; end

        def sin;   unop(Math.method(:sin)); end
        def cos;   unop(Math.method(:cos)); end
        def tan;   unop(Math.method(:tan)); end
        def asin;  unop(Math.method(:asin)); end
        def acos;  unop(Math.method(:acos)); end
        def atan;  unop(Math.method(:atan)); end
        def sinh;  unop(Math.method(:sinh)); end
        def cosh;  unop(Math.method(:cosh)); end
        def tanh;  unop(Math.method(:tanh)); end
        def asinh; unop(Math.method(:asinh)); end
        def acosh; unop(Math.method(:acosh)); end
        def atanh; unop(Math.method(:atanh)); end

        def sin!;   unop!(Math.method(:sin)); end
        def cos!;   unop!(Math.method(:cos)); end
        def tan!;   unop!(Math.method(:tan)); end
        def asin!;  unop!(Math.method(:asin)); end
        def acos!;  unop!(Math.method(:acos)); end
        def atan!;  unop!(Math.method(:atan)); end
        def sinh!;  unop!(Math.method(:sinh)); end
        def cosh!;  unop!(Math.method(:cosh)); end
        def tanh!;  unop!(Math.method(:tanh)); end
        def asinh!; unop!(Math.method(:asinh)); end
        def acosh!; unop!(Math.method(:acosh)); end
        def atanh!; unop!(Math.method(:atanh)); end


        def sqrt!; RAtlas::sqrt!(@storage); return self; end
        def sqrt; return clone.sqrt!; end
        def log!; RAtlas::log!(@storage); return self; end
        def log; return clone.log!; end
        def log10!; unop!(Math.method(:log10)); end
        def log10;  unop(Math.method(:log10)); end
        def exp!;   unop!(Math.method(:exp)); end
        def exp;    unop(Math.method(:exp)); end
        def abs!; map!{|f| f.abs}; end
        def abs; map{|f| f.abs}; end

        def coerce(arg)
            case arg
            when Numeric
                [Scalar.new(arg), self]
            when Scalar
                [arg.value, self]
            else
                raise TypeError, 
                "#{self.class} can't be coerced into #{other.class}"
            end
        end

        #
        # Sets elements in matrix according to position.  Numbering is
        # done row-wise.
        # 
        def indgen!(start = 0.0)
            RAtlas::indgen!(@storage, start)
            return self
        end
        
        # 
        # Copy the object and its vector/matrix.
        # 
        def clone
            return self.class.new(RAtlas::clone(@storage))
        end
        
        #
        # Copy the object, but dont care about vector/matrix elements.
        # 
        def dup; return self.class.new(RAtlas::dup(@storage)); end
        def to_s
            arr = to_a
            str = "#{self.class}["
            arr[0..-2].each {|row| 
                str += "[" + row.join(", ") + "], "
            }
            str += "[" + arr[-1].join(", ") + "]]"
            return str
        end
        def inspect; to_s; end    
        def size; RAtlas::size(@storage); end
        def column_size; size[1]; end
        def row_size; size[0]; end
        def get(idx)
            case idx
            when Range
                self.class.new(RAtlas::get_by_range(@storage,
                                       idx))
            when Array
                self.class.new(RAtlas::get_many(@storage, idx))
            else
                RAtlas::get_one(@storage, idx)
            end
        end

        #
        # Return the Fixnum value at index (i,j) for a Matrix.
        #
        def get2d_one(row, col)
            RAtlas::get2d_one(@storage, row, col)
        end

        #
        # Return the Matrix for the indexes (i,j) where i and j  an be
        # Fixnum, Array or Range.
        #
        def get2d_many(rows, cols)
            case rows
            when Fixnum
                case cols
                when Fixnum
                    RAtlas::get2d_one(@storage, rows, cols)
                when Array
                    self.class.new RAtlas::get2d_many(
                          @storage, [rows], cols)
                when Range
                    self.class.new RAtlas::get2d_by_range(
                     @storage, Range.new(rows, rows), cols)
                end
            when Range
                case cols
                when Range
                    self.class.new RAtlas::get2d_by_range(
                          @storage, rows, cols)
                when Array
                    rows = range2arr(rows, size[0])
                    self.class.new RAtlas::get2d_many(
                              @storage, rows, cols)
                when Fixnum
                    self.class.new RAtlas::get2d_by_range(
                      @storage, rows, Range.new(cols, cols))
                end
            when Array
                case cols
                when Range
                    cols = range2arr(cols, size[1])
                    self.class.new RAtlas::get2d_many(
                              @storage, rows, cols)
                when Array
                    self.class.new RAtlas::get2d_many(
                          @storage, rows, cols)
                when Fixnum
                    self.class.new RAtlas::get2d_many(
                        @storage, rows, [cols])
                end
            end

        end

        #
        # Set the value at index(es) idx to value val.  Index can be 
        # Fixnum, Array or Range.
        #
        def set!(idx, val)
            case idx
            when Range
                RAtlas::set_by_range!(@storage, idx, val.storage)
            when Array
                RAtlas::set_many!(@storage, idx, val.storage)
            else 
                RAtlas::set_one!(@storage, idx, val)
            end
            return self
        end
        def set2d!(rows, cols, vals)
            case vals
            when Fixnum
                set2d_one!(rows, cols, vals)
            else
                set2d_many!(rows, cols, vals)
            end
        end
        def set2d_one!(row, col, val)
            RAtlas::set2d_one!(@storage, row, col, val)
            return self
        end
        def set2d_many!(rows, cols, vals)
            case rows
            when Fixnum
                case cols
                when Fixnum
                    RAtlas::set2d_one!(@storage, rows, cols,
                              vals)
                when Array
                    RAtlas::set2d_many!(@storage, [rows],
                               cols, vals.storage)
                when Range
                    RAtlas::set2d_by_range!(@storage,
                        Range.new(rows, rows), cols,
                         vals.storage)
                end
            when Range
                case cols
                when Range
                    RAtlas::set2d_by_range!(@storage, rows,
                        cols, vals.storage)
                when Array
                    rows = range2arr(rows, size[0])
                    RAtlas::set2d_many!(@storage, rows, cols,
                               vals.storage)
                when Fixnum
                    RAtlas::set2d_by_range!(@storage, rows,
                      Range.new(cols, cols), vals.storage)
                end
            when Array
                case cols
                when Range
                    cols = range2arr(cols, size[1])
                    RAtlas::set2d_many!(@storage, rows, cols,
                               vals.storage)
                when Array
                    RAtlas::set2d_many!(@storage, rows, cols,
                           vals.storage)
                when Fixnum
                    RAtlas::set2d_many!(@storage, rows,
                               [cols], vals.storage)
                end
            end
            return self
        end
        def map!(&block)
            RAtlas::map!(@storage, block)
            return self
        end
        def map(&block); return clone.map!(&block); end
        def zip!(*args, &block)
            arg = args.map{|i| i.storage}
            RAtlas::zip!(@storage, arg, block)
            return self
        end
        def zip(*args, &block); return clone.zip!(*args, &block); end
        def col(idx)
            case idx
            when Fixnum
                Matrix.new(RAtlas::column(@storage, idx))
            when Array
                Matrix.new(RAtlas::columns(@storage, idx))
            else
                raise ArgumentError
            end
        end
        def cols(idx)
            Matrix.new(RAtlas::columns(@storage, idx))
        end
        def row(idx)
            case idx
            when Fixnum
                Matrix.new(RAtlas::row(@storage, idx))
            when Array
                Matrix.new(RAtlas::rows(@storage, idx))
            else
                raise ArgumentError
            end
        end
        def rows(idx); Matrix.new(RAtlas::rows(@storage, idx)); end
        def set_col!(idx, val)
            RAtlas::set_column(@storage,idx, val)
            return self
        end
        def set_cols!(idx, val)
            RAtlas::set_columns!(@storage, idx, val)
            return self
        end
        def set_row!(idx, val)
            RAtlas::set_row!(@storage, idx, val)
            return self
        end
        def set_rows!(idx, val)
            RAtlas::set_rows!(@storage, idx, val)
            return self
        end
        def re!; RAtlas::re!(@storage); return self; end
        def im!; RAtlas::im!(@storage); return self; end
        def re; self.class::new(RAtlas::re(@storage)); end
        def im; self.class::new(RAtlas::im(@storage)); end
        def sum; RAtlas::sum(@storage); end
        def rem!(arg); RAtlas::rem!(@storage, arg); return self; end
        def rem(arg); clone.rem!(arg); end
        def dotop!(scalop, matrop, arg, t = RAtlas::NOTRANS, 
               targ = RAtlas::NOTRANS)
            case arg
            when Numeric then scalop.call(@storage, arg)
            when RNum::Base
                matrop.call(t, @storage, targ, arg.storage)
            end
            return self
        end
        def dotop(scalop, matrop, arg, t = RAtlas::NOTRANS, 
              targ = RAtlas::NOTRANS)
            case arg
            when Numeric 
                self.class.new(scalop.call(@storage, arg))
            when RNum::Base
                self.class.new(matrop.call(t, @storage, targ, 
                               arg.storage))
            end
        end    
        def add!(*args)
            case args[0]
            when Numeric, self.class
                dotop!(RAtlas::method(:add!),
                       RAtlas::method(:madd!), *args)
            else
                raise ArgumentError, 
                    "Can only add #{self.class} or Numeric" 
                   end
        end
        def add(*args)
            case args[0]
            when Numeric, self.class
                dotop(RAtlas::method(:add), 
                      RAtlas::method(:madd), *args)
            else
                raise ArgumentError, 
                    "Can only add #{self.class} or Numeric" 
                   end
        end
        def sub!(*args)
            case args[0]
            when Numeric, self.class
                dotop!(RAtlas::method(:sub!),
                       RAtlas::method(:msub!), *args)
            else 
                raise ArgumentError,
                "Can only subtract #{self.class} or Numeric."
            end
        end
        def sub(*args)
            case args[0]
            when Numeric, self.class
                dotop(RAtlas::method(:sub),
                      RAtlas::method(:msub), *args)
            else 
                raise ArgumentError,
                "Can only subtract #{self.class} or Numeric."
            end
        end

        #
        # Scalar or element by element multiplication.  Equivalent with
        # .* found in numerical software.
        #
        def mul!(*args)
            dotop!(RAtlas::method(:mul!), RAtlas::method(:mmul!),
                   *args)
        end
        def mul(*args)
            dotop(RAtlas::method(:mul), RAtlas::method(:mmul),
                  *args)
        end

        #
        # Scalar or element by element division.  Equivalent with ./
        # found in numerical software.
        #
        def div!(*args)
            dotop!(RAtlas::method(:div!), RAtlas::method(:mdiv!),
                   *args)
        end
        def div(*args)
            dotop(RAtlas::method(:div), RAtlas::method(:mdiv),
                  *args)
        end

        def pow!(arg)
            RAtlas::pow!(@storage, arg)
            return self
        end
        def pow(arg); return clone.pow!(arg);end
        def **(arg); return clone.pow!(arg); end
        
        #
        # Returns sum_i sum_j op(A)ij * Bij
        #
        def colon(arg, t = RAtlas::NOTRANS, targ = RAtlas::NOTRANS)
            return RAtlas::colon(t, @storage, targ, arg.storage)
        end
        def +(arg); add(arg); end;
        def -(arg); sub(arg); end;
        def -@
            self.clone.mul!(-1.0)
        end
        def +@
            return self
        end
        def p(format = '%12.4e')
            puts pretty(format)
        end
        def pretty_print(arg)
            print pretty('%12.4e')
        end
        def concat(arg)
            Vector.new(RAtlas::concat(@storage, arg.storage))
        end
        def hcat(arg)
            Matrix.new(RAtlas::hcat(@storage, arg.storage))
        end
        def vcat(arg)
            Matrix.new(RAtlas::vcat(@storage, arg.storage))
        end
        def range2arr(range, len)
            fi = range.first
            la = range.last
            if la < 0
                   la += len
            end
             if range.exclude_end?
                la -= 1
            end
            Range.new(fi, la).to_a        
        end    
        
        #
        # Solve system AX = B.  arg <- X, A <- lu
        # 
        def solve!(arg)
            info, lustore, piv, xstore = 
                Lapack::gesv!(arg.storage, @storage)
            if info > 0
                sp = info - 1
                raise RNum::SingularError, \
                "Element (#{sp}, #{sp}) is exactly singular."
            end
            return self
        end
        def solve(arg); clone.solve!(arg.clone); end

        #
        # Back-substitution given a LU-matrix.  See +gertf+.
        #
        def getrs!(a, piv, t = RAtlas::NOTRANS)
            info, stor = Lapack::getrs!(t, a.storage, piv, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def getrs(*args); return clone.getrs!(*args); end

        #
        # Same as getrs, but matrix is symmetric positive definite.
        #
        def potrs!(a, uplo = RAtlas::UPPER)
            info, stor = Lapack::potrs!(uplo, a.storage, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def potrs(*args); clone.potrs(*args); end

        #
        # Special symmetric matrix back substitution.  See sytrf.
        #
        def sytrs!(a, piv, uplo = RAtlas::UPPER)
            info, stor = Lapack::sytrs!(uplo, a.storage, piv,
                            @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def sytrs(*args); return clone.sytrs!(*args); end
    
        #
        # Inner product or scalar multiplication.
        #
        def *(arg)
            case arg
            when Numeric
                mul(arg)
            when RNum::Base
                inner(arg)
            else
                raise ArgumentError, "Can not multiply #{self.class} with #{arg.class}"
            end
        end
        def <(arg)
          RNum::compare(self, :<, arg)
        end
        def >(arg)
          RNum::compare(self, :>, arg)
        end
        def ==(arg)
          RNum::compare(self, :==, arg)
        end
        def <=(arg)
          RNum::compare(self, :<=, arg)
        end
        def >=(arg)
          RNum::compare(self, :>=, arg)
        end
    end


    class Matrix < Base
        public_class_method :new
        def Matrix::[](*args)
            case args[0]
            when RNum::Base then new(args[0])
            when Array then new(args)
            else
                raise ArgumentError,
                "Expected Matrix, Vector or Numeric"
            end
        end
        def Matrix::zeros(arg)
            new(RAtlas::zeros(arg))
        end
        def Matrix::ones(arg)
            new(RAtlas::ones(arg))
        end
        def Matrix::eye(arg)
            new(RAtlas::eye(arg))
        end
        def Matrix::alloc(arg)
            new(arg)
        end
        def initialize(arg)
            case arg
                when RAtlas::Storage then @storage = arg
                when Array then @storage = RAtlas::matrix(arg)
                when RNum::Base
                    @storage = arg.clone.storage
                else raise ArgumentError
            end
            return self
        end
        def to_a; RAtlas::matrix_to_a(@storage); end
        def to_v 
            if size[1] == 1
                Vector::new(RAtlas::clone(@storage))
            else
                raise ArgumentError, 
                "Can only cast 1-column matrix to Vector."
            end
        end
        # 
        # Matrix pretty-print.
        #
        def pretty(format = '%12.4e')
            str = String.new()
            to_a.each{|a| 
                a.each{|ai| 
                    str += sprintf(format, ai)
                }
                str += "\n"
            }
            return str
        end
        #
        # Return a string representation of the matrix with insidences of
        # negative, positive and zeros elements as {I, -I, -, +, .}.
        def spy(tol=Float::EPSILON)
            str = String.new()
            self.to_a.each{|row|
                str += row.map{|e| 
                    case e
                    when -1.0 then 'N'
                    when 1.0 then 'I'
                    else
                        if e.abs < tol
                            '.'
                        else
                            e < 0.0 ? '-' : '+'
                        end
                    end
                }.to_s + "\n"
            }
            return str
        end
        #
        # Change shape to arg dimensions.  Elements are read
        # column-wise.  RNum::eye([2,3]).reshape!([3,2])
        #   => RNum::Matrix[[1.0, 1.0], [0.0, 0.0], [0.0, 0.0]]
        #   
        def reshape!(arg)
            RAtlas::reshape!(@storage, arg)
            return self
        end
        def reshape(arg); clone.reshape!(arg); end
        def flatten!
            nelem = size[0]*size[1]
            RAtlas::reshape!(@storage, [1, nelem])
            return self
        end
        def flatten; clone.flatten!;  end
        def [](r, c); get2d_many(r, c); end
        def []=(r, c, val); set2d_many!(r, c, val); end
        
        #
        # Returns matrix transpose.
        # 
        def t; self.class.new(RAtlas::transpose(@storage)); end
        def t! 
            @storage = RAtlas::transpose(@storage)
            return self
        end
        def d; Vector.new(RAtlas.diag2vec(@storage)); end
        def /(arg)
            case arg
            when Numeric
                div(arg)
            when RNum::Matrix
                solve(arg)
            else
                raise ArgumentError, 
                    "/-operator for #{self.class} only."
            end
        end
        def dmul!(arg)
            case arg
            when RNum::Vector
                RAtlas::admul!(@storage, arg.storage)
                return self
            else
                raise ArgumentError, "Expected Vector."
            end
        end
        def dmul(arg); clone.dmul!(arg); end
        def dadd!(arg)
            case arg
            when RNum::Vector
                RAtlas::dadd!(@storage, arg.storage)
                return self
            else
                raise ArgumentError, "Expected Vector."
            end
        end
        def dadd(arg); clone.dadd!(arg); end
        def dsub!(arg)
            case arg
            when RNum::Vector
                RAtlas::dsub!(@storage, arg.storage)
                return self
            else
                raise ArgumentError, "Expected Vector."
            end
        end
        def dsub(arg); clone.dsub!(arg); end
        
        #
        # Sum over rows.
        #
        def rowsum; self.class.new(RAtlas::rowsum(@storage)); end
        
        #
        # Sum over columns.
        #
        def colsum; self.class.new(RAtlas::colsum(@storage)); end

        #
        # Compute the singular values of a general M-by-N matrix A
        #
        def singular_values
            info, astor, sstor, ustor, vtstor =
                Lapack::gesvd!(clone.storage)
            if info > 0
                raise RNum::Error, "Algorithm did not converge."
            end
            RNum::Vector.new sstor
        end
        
        # 
        # Compute the rank by counting the number of positive
        # singular values, and taking into account they are
        # ordered non-increasingly.
        #
        def rank(tol=Float::EPSILON*10e7)
            (singular_values.to_a << 0).each_with_index \
            do |v,i|
                return i if v.abs < tol
            end
        end

        #
        # Computes an LU factorization of a general M-by-N matrix A
        # using partial pivoting with row interchanges.
        #
        # The factorization has the form
        #    A = P * L * U
        # where P is a permutation matrix, L is lower triangular with
        # unit diagonal elements (lower trapezoidal if m > n), and U is
        # upper triangular (upper trapezoidal if m < n).
        # 
        # Row i of the matrix was interchanged with row piv[i].
        # If info > 0, then the element U(i,i) is exactly singular.
        # 
        def getrf!
            info, stor, piv = Lapack::getrf!(@storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return [self, piv]
        end
        def getrf; return clone.getrf!; end

        #
        # Same as lu, but matrix is symmetric positive definite.
        #
        def potrf!(uplo = RAtlas::UPPER)
            info, stor = Lapack::potrf!(uplo, @storage)
            if info > 0
                raise RNum::Error, \
                "The leading minor of order #{info} is not positive definite."
            end
            return self
        end
        def potrf(*args); clone.chol!(*args); end
        alias :chol! potrf!
        alias :chol potrf

        #
        # Factorization for generic symmetrical matrix. A is special
        # matrix as described in Lapack function xsytrf (x = {d, z}).    
        #
        def sytrf!(uplo = RAtlas::UPPER)
            info, stor, piv = Lapack::sytrf!(uplo, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return [self, piv]
        end
        def sytrf(*args); clone.sytrf!(*args); end

        #
        # Matrix inverse.
        #
        def inv!
            info, stor, piv = Lapack::getrf!(@storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            info, stor = Lapack::getri!(@storage, piv)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def inv; return clone.inv!; end

        #
        # Same as inv, but matrix is symmetric positive definite.
        #
        def poinv!(uplo = RAtlas::UPPER)
            info, stor = Lapack::potrf!(uplo, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            info, stor = Lapack::potri!(uplo, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def poinv(*args); clone.poinv!(*args); end

        #
        # Same as inv, but matrix is symmetric.
        #
        def syinv!(uplo = RAtlas::UPPER)
            info, stor, piv = Lapack::sytrf!(uplo, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            info, stor = Lapack::sytri!(uplo, @storage, piv)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def syinv(*args); clone.syinv!(*args); end

        #
        # Solve system AX = B. arg <- X, A <- self.polu
        # A is a symmetric positive definite matrix.
        #
        def posolve!(arg, uplo = RAtlas::LOWER)
            info, lstore, xstore = 
                Lapack::posv!(uplo, arg.storage, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def posolve(arg, uplo = RAtlas::LOWER)
            clone.posolve!(arg.clone, uplo)
        end

        #
        # Solve system AX = B. arg <- X, A <- self.sylu
        # A is a symmetric matrix.
        #
        def sysolve!(arg, uplo = RAtlas::LOWER)
            info, lustore, piv, xstore = 
                Lapack::sysv!(uplo, arg.storage, @storage)
            if info > 0
                raise RNum::SingularError, \
                "Element (#{info-1}, #{info-1}) is exactly singular."
            end
            return self
        end
        def sysolve(arg, uplo = RAtlas::LOWER)
            clone.sysolve!(arg.clone, uplo)
        end

        #
        # Inner product <self, arg>
        #
        def inner(arg, t = RAtlas::NOTRANS, targ = RAtlas::NOTRANS)
            case arg
            when RNum::Vector
                t == RAtlas::TRANS ?  m = size[1] : m = size[0]
                c = RNum::alloc(m)
                Blas::gemv!(t, 1.0, @storage, arg.storage, 0.0,
                    c.storage)
                return c
            when RNum::Matrix
                t == RAtlas::TRANS ? m = size[1] : m = size[0]
                if targ == RAtlas::TRANS
                    n = arg.size[0]
                else
                    n = arg.size[1]
                end
                c = RNum::alloc([m,n])
                Blas::gemm!(t, targ, 1.0, @storage, 
                            arg.storage, 0.0, c.storage)
                return c
            end
        end

        #
        # Dyad or outer product.
        # 
        def dyad(arg, t = RAtlas::NOTRANS, targ = RAtlas::NOTRANS)
            case arg
            when RNum::Matrix
                t == RAtlas::TRANS ? m = size[1] : m = size[0]
                if targ == RAtlas::TRANS
                    n = arg.size[1]
                    targ = RAtlas::NOTRANS
                else
                    n = arg.size[0]
                    targ = RAtlas::TRANS
                end
                c = RNum::alloc([m,n])
                Blas::gemm!(t, targ, 1.0, @storage,
                        arg.storage, 0.0, c.storage)
                return c
            else
                raise ArgumentError, "Expect matrix."
            end
        end
        def dyad!(*args)
            c = dyad(*args)
            @storage = c.storage
            return self
        end

        #
        # Dyad or outer product (A*A' in numerical software).
        #
        def dyadic(t = RAtlas::NOTRANS)
            t == RAtlas::TRANS ? m = size[1] : m = size[0]
            c = RNum::alloc([m,m])
            Blas::syrk!(t, 1.0, @storage, 0.0, c.storage)
            return c
        end
        def dyadic!(*args)
            c = dyadic(*args)
            @storage = c.storage
            return self
        end

        # 
        # Rank2 dyad (A*B' + B*A' in numerical software).
        # 
        def dyadic2(arg, t = RAtlas::NOTRANS, targ = RAtlas::NOTRANS) 
            if t == RAtlas::TRANS
                astor = RAtlas::transpose(@storage)
            else
                astor = @storage
            end
            if targ == RAtlas::TRANS
                bstor = RAtlas::transpose(arg.storage)
            else
                bstor = arg.storage
            end
            c = RNum::alloc([size[0], arg.size[0]])
            Blas::syr2k!(RAtlas::NOTRANS, 1.0, astor, bstor, 0.0,
                     c.storage)
            return c
        end
        def dyadic2!(*args)
               c = dyadic2(*args)
               @storage = c.storage
               return self
        end
        def to_gplot (x = nil, y = nil)
            xgrid = x || (0...self.size[1]).to_a
            ygrid = y || (0...self.size[0]).to_a
            f = ""
            ygrid.length.times do |j|
                y = ygrid[j]
                xgrid.length.times do |i|
                    if ( self[j,i] )
                        f << "#{xgrid[i]} #{y} #{self[j,i]}\n"
                    end
                end
            end
            f
        end
        def to_gsplot; to_gplot; end
    end




    class Vector < Base
        public_class_method :new
        def Vector::[](*args)
            case args[0]
            when RNum::Base then new(args[0])
            when Numeric then new(args)
            when Array then new(*args)
            else
                raise ArgumentError, "Expected Matrix or Numeric"
            end
        end
        def Vector::zeros(arg)
            return new(RAtlas::zeros([arg, 1]))
        end
        def Vector::ones(arg)
            return new(RAtlas::ones([arg,1]))
        end
        def Vector::eye(arg)
            return new(RAtlas::eye([arg,1]))
        end
        def Vector::alloc(arg)
            new(arg)
        end

        def initialize(arg)
            case arg
            when RAtlas::Storage then @storage = arg
            when Array then @storage = RAtlas::vector(arg)
            when RNum::Vector
                @storage = RAtlas::clone(arg.storage)
            when RNum::Base
                case arg.size[-1]
                when 1
                    @storage = RAtlas::clone(arg.storage)
                else
                    @storage = arg.col(0).storage
                end
            else
                raise ArgumentError
            end
            return self
        end

        #
        # Cast vector to matrix with vector elements on diagonal.
        #
        def d; Matrix.new(RAtlas.vec2diag(@storage)); end

        #
        # Cast vector to 1-row matrix.
        #
        def row; reshape(1, length); end

        #
        # Cast vector to 1-column matrix.
        #
        def col; reshape(length, 1); end
        def reshape(*arg)
            case arg.size
            when 1
                Matrix.new(RAtlas::reshape!(clone.storage, 
                                arg[0]))
            when 2
                Matrix.new(RAtlas::reshape!(clone.storage,
                                arg))
            else
                raise ArgumentError
            end
        end
        def []=(*args)
            val = args.last
            idx = args.length > 2 ? args[0..-2] : args[0]
            set!(idx, val)
        end
        def [](*i); i.size == 1 ? get(i[0]) : get(i); end

        #
        # self.dmul!(a): Eqivalent to diag(v)*A in numerical software.
        #
        def dmul!(arg)
            case arg
            when RNum::Matrix
                RAtlas::damul!(arg.storage, @storage)
                return arg
            else
                raise ArgumentError, "Expected Matrix."
            end
        end
        def dmul(arg); dmul!(arg.clone); end
        def length; RAtlas::size(@storage)[0]; end
        def to_a; RAtlas::vector_to_a(@storage); end
        def to_s
            arr = to_a
            str = "#{self.class}["
            arr[0..-2].each {|elem| 
                str += elem.to_s + ", " 
            }
            str += arr[-1].to_s + "]"
            return str
        end

        #
        # Return string for pretty-printing.
        #
        def pretty(format = '%12.4e')
            str = String.new()
            to_a.each{|i| str += sprintf(format, i)} 
            return str
        end
        # 
        # Transpose
        #
        def t; self; end
        def t!; self; end

        #
        # Sum of absolute value for all elements.
        # 
        def asum; Blas::asum(@storage); end

        #
        # Index of absolute maximum value.
        #
        def iamax; Blas::iamax(@storage); end

        #
        # Maximum absolute value.
        #
        def amax; get(iamax).abs; end
    
        #
        # Minimum absoulte value
        #
        def amin; abs.min; end

        #
        # Maximum value.
        #
        def max(arg=nil)
               case arg
               when nil
                   to_a.max
               when Numeric
                   [to_a.max, arg.to_f].max
               else
                   raise ArgumentError, "Expect numeric."
               end
        end

        #
        # Minimum value.
        #
        def min(arg=nil)
            case arg
            when nil
                to_a.min
            when Numeric
                [to_a.min, arg.to_f].min
            else
                   raise ArgumentError, "Expect numeric."
            end
        end

        # Norm.  Default is 2-norm.
        def norm(p=2.0)
            if p.to_f == 2.0
                Blas::nrm2(@storage)
            else
                case p.infinite?
                when 1
                    self[amax]
                when -1
                    abs.min
                else
                    map{|i| i**p}.asum**(1.0/p)
                end
            end
        end

        #
        # Scalar division.
        #
        def /(arg)
            case arg
            when Numeric
                div(arg)
            else
                raise ArgumentError, 
                "/-operator on #{self.class} for Numeric only."
            end
        end

        #
        # Inner product <self, arg>
        #
        def inner(arg, t = RAtlas::NOTRANS) # a'*B
            case arg
            when RNum::Vector
                Blas::dot(@storage, arg.storage)
            when RNum::Matrix
            #    Blas has <M,v>. Get <v, M> from <a, B>=<B',a>'
                if t == RAtlas::TRANS
                    t = RAtlas::NOTRANS
                else
                    t = RAtlas::TRANS
                end
                n = self.length
                retval = RNum::alloc(n)
                Blas::gemv!(t, 1.0, arg.storage, @storage, 0.0, 
                                              retval.storage)
                n == 1 ? retval.get(0) : retval
            else
                raise ArgumentError, "Unable to calculate inner product of #{self.class} with #{arg.class}"
            end
        end

        #
        # Dyadic or outer product A = x*y'
        #
        def dyad(arg)
            case arg
            when RNum::Vector
                n = arg.length
                a = RNum.zeros([n,n])
                Blas::ger!(1.0, @storage, arg.storage,
                       a.storage)
                return a
            else
                raise ArgumentError,
                 "Dyad for vector with vector only."
            end
        end

        #
        # Dyadic or outer product A = x*x'
        #
        def dyadic
            n = length
            a = RNum.zeros([n,n])
            Blas::syr!(1.0, @storage, a.storage)
            return a
        end
        
        #
        # Dyadic rank2 outer product A = self*arg' + arg*self'
        #
        def dyadic2(arg)
            case arg
            when RNum::Vector
                n = arg.length
                a = RNum.zeros([n,n])
                Blas::syr2!(1.0, @storage, arg.storage,
                        a.storage)
                return a
            else
                raise ArgError, "Dyadic2 for vector only."
            end
        end
        def to_gplot
            s = ""
            self.length.times{ |i| s << "#{self[i]}\n" }
            s
        end

    end




    class Scalar < Numeric # :nodoc:
        # Need to implement **-function for matrix exponentials.
        attr_reader :value
        def initialize(value)
            @value = value
        end
        def +(other)
            case other
            when Numeric
                @value + other
            when RNum::Scalar
                @value + other.value
            when RNum::Base
                other.add(@value)
            else
                x, y = other.coerce(@value)
                x + y
            end
        end
        def -(other)
            case other
            when Numeric
                @value - other
            when RNum::Scalar
                @value - other.value
            when RNum::Base
                (-other).add(@value)
            else
                x, y = other.coerce(@value)
                x + y
            end
        end
        def *(other)
            case other
            when Numeric
                @value * other
            when RNum::Scalar
                @value  * other.value
            when RNum::Base
                other.mul(@value)
            else
                x, y = other.coerce(@value)
                x * y
            end
        end
        def /(other)
            case other
            when Numeric
                @value / other
            when RNum::Scalar
                @value / other.value
            when RNum::Vector
                (other**-1.0).mul(@value)
            when RNum::Base
                other.inv.mul(@value)
            else
                x, y = other.coerce(@value)
                x / y
            end
        end
        def coerce(other)
            case other
            when Numeric
                [other, @value]
            else
                raise TypeError, 
                "#{self.class} can't be coerced into #{other.class}"
            end
        end
    end


    class Error < StandardError; end
    class SingularError < RNum::Error; end

    # to be implemented:
    #class RNum::Diagonal < RNum::Matrix; end
    #class RNum::Cmatrix < RNum::Matrix; end
    #class RNum::CVector < RNum::Vector; end
end



