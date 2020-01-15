package fields;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import fields.PrimeField.PrimeFieldElement;

public class PrimeField extends AbstractField<PrimeFieldElement> {
	private int p;
	private Map<Integer,PrimeFieldElement> elements;

	public PrimeField(int p) {
		this.p = p;
		this.elements = new TreeMap<Integer, PrimeFieldElement>();
		if (!BigInteger.valueOf(p).isProbablePrime(10))
			throw new RuntimeException("Not a prime number");
	}

	@Override
	public PrimeFieldElement zero() {
		return this.getElement(0);
	}

	@Override
	public PrimeFieldElement one() {
		return this.getElement(1);
	}

	@Override
	public int characteristic() {
		return p;
	}

	@Override
	public PrimeFieldElement add(PrimeFieldElement t1, PrimeFieldElement t2) {
		return this.getElement(t1.value + t2.value);
	}

	@Override
	public PrimeFieldElement negative(PrimeFieldElement t) {
		return this.getElement(-t.value);
	}

	@Override
	public PrimeFieldElement multiply(PrimeFieldElement t1, PrimeFieldElement t2) {
		return this.getElement(t1.value * t2.value);
	}

	@Override
	public PrimeFieldElement inverse(PrimeFieldElement t) {
		int u, v, x, q, r, y;
		u = this.p;
		v = t.value;
		y = 0;
		x = 1;
		while (u % v != 0) {
			q = u / v;
			r = u % v;
			u = v;
			v = r;
			r = x;
			x = y - q * x;
			y = r;
		}
		return this.getElement(x);
	}

        private int jacobiSymbol(int a, int b) {
          if (b % 2 == 0) {
            throw new RuntimeException("JacobiSymbol (" + a + "/" + b + ") undefined");
          }
          a = a % b;
          if (a + 1 == b) {
            if (b % 4 == 1) {
              return 1;
            } else if (b % 4 == 3) {
              return -1;
            }
          }
          if (a == 2) {
            if ((b % 8 == 1) || (b % 8 == 7))
              return 1;
            else if ((b % 8 == 3) || (b % 8 == 5)) {
              return -1;
            }
          }
          if (a % 2 == 0) {
            return jacobiSymbol(2, b) * jacobiSymbol(a / 2, b);
          }
          int inv = jacobiSymbol(b, a);
          if (a % 4 == 3 && b % 4 == 3) {
            return -inv;
          }
          return inv;
        }

        @Override
        public boolean hasSqrt(PrimeFieldElement t) {
          if (this.p == 2) {
            return true;
          }
          return jacobiSymbol(t, this.p) == 1;
        }

        @Override
        public PrimeFieldElement sqrt(PrimeFieldElement t) {
          if (this.p % 4 == 3) {
            return this.power(t, (this.p + 1) / 4);
          }
          if (this.p == 2) {
            return t;
          }
          if (t.equals(this.zero())) {
            return this.zero();
          }
          int s = 0;
          int q = this.p - 1;
          while (q % 2 == 0) {
            s++;
            q /= 2;
          }
          int z = 2;
          while (jacobiSymbol(z, this.p) == 1) {
            z++;
          }
          int m = s;
          PrimeFieldElement c = this.power(z, q);
          PrimeFieldElement n = this.power(t, q);
          PrimeFieldElement r = this.power(t, (q + 1) / 2);
          while (!n.equals(this.one())) {
            int i = 0;
            PrimeFieldElement sq = n;
            while (!sq.equals(this.one())) {
              sq = this.multiply(sq, sq);
              i++;
            }
            PrimeFieldElement b = c;
            for (int j = 0; j < m - i - 1) {
              b = this.multiply(b, b);
            }
            m = i;
            c = this.multiply(b, b);
            n = this.multiply(n, c);
            r = this.multiply(r, b);
          }
          return r;
        }

	@Override
	public PrimeFieldElement getRandomElement() {
		return this.getElement((int)Math.floor((Math.random() * this.p)));
	}

	@Override
	public int getNumberOfElements() {
		return this.p;
	}

	@Override
	public Iterable<PrimeFieldElement> getElements() throws InfinityException {
		return new Iterable<PrimeFieldElement>() {
			
			@Override
			public Iterator<PrimeFieldElement> iterator() {
				return new Iterator<PrimeFieldElement>() {
					private int i = 0;
					@Override
					public boolean hasNext() {
						return this.i < p;
					}

					@Override
					public PrimeFieldElement next() {
						this.i++;
						return PrimeField.this.getElement(i - 1);
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};
	}
	@Override
	public boolean isFinite() {
		return true;
	}
	@Override
	public String toString() {
		return "F" + p;
	}
	public PrimeFieldElement getElement(int value) {
		if (this.elements.containsKey(value))
			return this.elements.get(value);
		this.elements.put(value, new PrimeFieldElement(value));
		return this.getElement(value);
	}
	public class PrimeFieldElement implements Element {
		private int value;
		
		private PrimeFieldElement(int value) {
			this.value = value % p;
			while (this.value < 0)
				this.value += p;
		}

		public int getValue() {
			return value;
		}
		
		public boolean equals(Object o) {
			if (!(o instanceof PrimeFieldElement))
				return false;
			PrimeFieldElement e = (PrimeFieldElement)o;
			return this.value == e.value;
		}
		public String toString() {
			return Integer.toString(this.value);
		}

		@Override
		public int compareTo(Element e) {
			PrimeFieldElement o = (PrimeFieldElement)e;
			return this.value - o.value;
		}
	}
}
