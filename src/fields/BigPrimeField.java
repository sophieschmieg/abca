package fields;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.Random;

import fields.BigPrimeField.BigPrimeFieldElement;

public class BigPrimeField extends AbstractField<BigPrimeFieldElement> {
	private BigInteger p;

	public BigPrimeField(BigInteger p) {
		this.p = p;
		if (!p.isProbablePrime(10))
			throw new RuntimeException("Not a prime number");
	}

	@Override
	public BigPrimeFieldElement zero() {
		return this.getElement(BigInteger.ZERO);
	}

	@Override
	public BigPrimeFieldElement one() {
		return this.getElement(BigInteger.ONE);
	}

	@Override
	public int characteristic() {
		return p.intValue();
	}

	@Override
	public BigPrimeFieldElement add(BigPrimeFieldElement t1, BigPrimeFieldElement t2) {
		return this.getElement(t1.value.add(t2.value));
	}

	@Override
	public BigPrimeFieldElement negative(BigPrimeFieldElement t) {
		return this.getElement(t.value.negate());
	}

	@Override
	public BigPrimeFieldElement multiply(BigPrimeFieldElement t1, BigPrimeFieldElement t2) {
		return this.getElement(t1.value.multiply(t2.value));
	}

	@Override
	public BigPrimeFieldElement inverse(BigPrimeFieldElement t) {
		BigInteger u, v, x, q, r, y;
		u = this.p;
		v = t.value;
		y = BigInteger.ZERO;
		x = BigInteger.ONE;
		while (!u.mod(v).equals(BigInteger.ZERO)) {
                        BigInteger[] qr = u.divideAndRemainder(v);
			q = qr[0];
			r = qr[1];
			u = v;
			v = r;
			r = x;
			x = y.subtract(q.multiply(x));
			y = r;
		}
		return this.getElement(x);
	}

	@Override
	public BigPrimeFieldElement getRandomElement() {
		return this.getElement(BigInteger.valueOf(new Random().nextInt()));
	}

	@Override
	public int getNumberOfElements() {
		return this.p.intValue();
	}

	@Override
	public Iterable<BigPrimeFieldElement> getElements() throws InfinityException {
		return new Iterable<BigPrimeFieldElement>() {
			
			@Override
			public Iterator<BigPrimeFieldElement> iterator() {
				return new Iterator<BigPrimeFieldElement>() {
					private BigInteger i = BigInteger.ZERO;
					@Override
					public boolean hasNext() {
						return this.i.compareTo(p) < 0;
					}

					@Override
					public BigPrimeFieldElement next() {
                                                BigInteger tmp = this.i;
						this.i = this.i.add(BigInteger.ONE);
						return BigPrimeField.this.getElement(tmp);
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
	public BigPrimeFieldElement getElement(BigInteger value) {
		return new BigPrimeFieldElement(value);
	}
	public class BigPrimeFieldElement implements Element {
		private BigInteger value;
		
		private BigPrimeFieldElement(BigInteger value) {
			this.value = value.mod(p);
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(p);
		}

		public BigInteger getValue() {
			return value;
		}
		
		public boolean equals(Object o) {
			if (!(o instanceof BigPrimeFieldElement))
				return false;
			BigPrimeFieldElement e = (BigPrimeFieldElement)o;
			return this.value.equals(e.value);
		}
		public String toString() {
			return this.value.toString();
		}

		@Override
		public int compareTo(Element e) {
			BigPrimeFieldElement o = (BigPrimeFieldElement)e;
			return this.value.compareTo(o.value);
		}
	}
}
