package fields;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.ModularIntegerRing.ModularIntegerRingElement;

public class ModularIntegerRing extends AbstractRing<ModularIntegerRingElement> {
	private int n;

	public ModularIntegerRing(int n) {
		this.n = n;
	}

	@Override
	public ModularIntegerRingElement zero() {
		return this.getElement(0);
	}

	@Override
	public ModularIntegerRingElement one() {
		return this.getElement(1);
	}

	@Override
	public int characteristic() {
		return n;
	}

	@Override
	public ModularIntegerRingElement add(ModularIntegerRingElement t1, ModularIntegerRingElement t2) {
		return this.getElement(t1.value + t2.value);
	}

	@Override
	public ModularIntegerRingElement negative(ModularIntegerRingElement t) {
		return this.getElement(-t.value);
	}

	@Override
	public ModularIntegerRingElement multiply(ModularIntegerRingElement t1, ModularIntegerRingElement t2) {
		return this.getElement(t1.value * t2.value);
	}

	@Override
	public ModularIntegerRingElement inverse(ModularIntegerRingElement t) {
		int u, v, x, q, r, y;
		u = this.n;
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

	@Override
	public ModularIntegerRingElement getRandomElement() {
		return this.getElement((int)Math.floor((Math.random() * this.n)));
	}

	@Override
	public int getNumberOfElements() {
		return this.n;
	}

	@Override
	public Iterable<ModularIntegerRingElement> getElements() throws InfinityException {
		return new Iterable<ModularIntegerRingElement>() {
			
			@Override
			public Iterator<ModularIntegerRingElement> iterator() {
				return new Iterator<ModularIntegerRingElement>() {
					private int i = 0;
					@Override
					public boolean hasNext() {
						return this.i < n;
					}

					@Override
					public ModularIntegerRingElement next() {
						this.i++;
						return ModularIntegerRing.this.getElement(i - 1);
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
		return "Z/" + n + "Z";
	}
	public ModularIntegerRingElement getElement(int value) {
		return new ModularIntegerRingElement(value);
	}
	public class ModularIntegerRingElement implements Element {
		private int value;
		
		private ModularIntegerRingElement(int value) {
			this.value = value % n;
			while (this.value < 0)
				this.value += n;
		}

		public int getValue() {
			return value;
		}
		
		public boolean equals(Object o) {
			if (!(o instanceof ModularIntegerRingElement))
				return false;
			ModularIntegerRingElement e = (ModularIntegerRingElement)o;
			return this.value == e.value;
		}
		public String toString() {
			return Integer.toString(this.value);
		}

		@Override
		public int compareTo(Element e) {
			ModularIntegerRingElement o = (ModularIntegerRingElement)e;
			return this.value - o.value;
		}
	}
	@Override
	public boolean isUnit(ModularIntegerRingElement t) {
		int u, v, r;
		u = this.n;
		v = t.value;
		while (u % v != 0) {
			r = u % v;
			u = v;
			v = r;
		}
		return v == 1;
	}

	@Override
	public int getNumberOfUnits() {
		Map<Integer,Integer> primes = new TreeMap<Integer,Integer>();
		int t = this.n;
		for (int i = 2; i <= t; i++) {
			if (t % i == 0) {
				int j = 0;
				while (t % i == 0) {
					j++;
					t /= i;
				}
				primes.put(i, j);
			}
		}
		int result = 1;
		for (int p : primes.keySet()) {
			for (int k = 0; k < primes.get(p); k++)
				result *= p;
			result *= p-1;
		}
		return result;
	}

	@Override
	public Iterable<ModularIntegerRingElement> getUnits() {
		return new Iterable<ModularIntegerRingElement>() {
			
			@Override
			public Iterator<ModularIntegerRingElement> iterator() {
				return new Iterator<ModularIntegerRingElement>() {
					private int i = 0;
					@Override
					public boolean hasNext() {
						return this.i < n;
					}

					@Override
					public ModularIntegerRingElement next() {
						do {
							this.i++;
						} while (!ModularIntegerRing.this.isUnit(ModularIntegerRing.this.getElement(i)));
						return ModularIntegerRing.this.getElement(i);
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
	public boolean isIntegral() {
		return BigInteger.valueOf(n).isProbablePrime(10);
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public List<ModularIntegerRingElement> quotientAndRemainder(ModularIntegerRingElement dividend, ModularIntegerRingElement divisor) {
		throw new RuntimeException();
	}
}
