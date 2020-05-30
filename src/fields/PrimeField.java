package fields;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import fields.PrimeField.PrimeFieldElement;
import util.MiscAlgorithms;

public class PrimeField extends AbstractField<PrimeFieldElement> {
	private BigInteger p;

	public PrimeField(BigInteger p) {
		this.p = p;
		if (!p.isProbablePrime(10))
			throw new RuntimeException("Not a prime number");
	}
	
	public PrimeField(int p) {
		this.p = BigInteger.valueOf(p);
		if (!this.p.isProbablePrime(10))
			throw new RuntimeException("Not a prime number");
	}

	@Override
	public PrimeFieldElement zero() {
		return this.getElement(BigInteger.ZERO);
	}

	@Override
	public PrimeFieldElement one() {
		return this.getElement(BigInteger.ONE);
	}

	@Override
	public BigInteger characteristic() {
		return p;
	}

	@Override
	public PrimeFieldElement add(PrimeFieldElement t1, PrimeFieldElement t2) {
		return this.getElement(t1.value.add(t2.value));
	}

	@Override
	public PrimeFieldElement negative(PrimeFieldElement t) {
		return this.getElement(t.value.negate());
	}

	@Override
	public PrimeFieldElement multiply(PrimeFieldElement t1, PrimeFieldElement t2) {
		return this.getElement(t1.value.multiply(t2.value));
	}

	@Override
	public PrimeFieldElement inverse(PrimeFieldElement t) {
		return this.getElement(t.value.modInverse(p));
	}

	private int jacobiSymbol(BigInteger a, BigInteger b) {
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger four = BigInteger.valueOf(4);
		BigInteger eight = BigInteger.valueOf(8);
		if (b.mod(two).equals(zero)) {
			throw new RuntimeException("JacobiSymbol (" + a + "/" + b + ") undefined");
		}
		a = a.mod(b);
		// a == -1 mod b
		int aModFour = a.mod(four).intValue();
		int bModFour = b.mod(four).intValue();
		if (a.add(one).equals(b)) {
			if (bModFour == 1) {
				return 1;
			} else if (bModFour == 3) {
				return -1;
			}
		}
		if (a.equals(two)) {
			int bModEight = b.mod(eight).intValue();
			if (bModEight == 1 || bModEight == 7)
				return 1;
			else if (bModEight == 3 || bModEight == 5) {
				return -1;
			}
		}
		if (aModFour % 2 == 0) {
			return jacobiSymbol(two, b) * jacobiSymbol(a.shiftRight(1), b);
		}
		int inv = jacobiSymbol(b, a);
		if (aModFour == 3 && bModFour == 3) {
			return -inv;
		}
		return inv;
	}

	@Override
	public boolean hasSqrt(PrimeFieldElement t) {
		if (this.p == BigInteger.TWO) {
			return true;
		}
		int jacobi = jacobiSymbol(t.value, this.p);
		return jacobi == 1 || jacobi == 0;
	}

	@Override
	public Set<PrimeFieldElement> sqrt(PrimeFieldElement t) {
		if (t.equals(this.zero())) {
			return Collections.singleton(this.zero());
		}
		if (jacobiSymbol(t.value, this.p) == -1) {
			return Collections.emptySet();
		}
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger four = BigInteger.valueOf(4);
		if (this.p.mod(four).intValue() == 3) {
			PrimeFieldElement sqrt = this.power(t, this.p.add(one).divide(four));
			Set<PrimeFieldElement> result = new TreeSet<>();
			result.add(sqrt);
			result.add(this.negative(sqrt));
			return result;
		}
		if (this.p.equals(two)) {
			return Collections.singleton(t);
		}
		int s = 0;
		BigInteger q = this.p.subtract(one);
		while (q.mod(two).equals(zero)) {
			s++;
			q = q.shiftRight(1);
		}
		BigInteger z = two;
		while (z.mod(this.p).equals(zero) || jacobiSymbol(z, this.p) == 1) {
			z = z.add(one);
		}
		int m = s;
		PrimeFieldElement c = this.power(this.getElement(z), q);
		PrimeFieldElement n = this.power(t, q);
		PrimeFieldElement r = this.power(t, q.add(one).shiftRight(1)); // (q + 1) / 2
		while (!n.equals(this.one())) {
			int i = 0;
			PrimeFieldElement sq = n;
			while (!sq.equals(this.one())) {
				sq = this.multiply(sq, sq);
				i++;
			}
			PrimeFieldElement b = c;
			for (int j = 0; j < m - i - 1; j++) {
				b = this.multiply(b, b);
			}
			m = i;
			c = this.multiply(b, b);
			n = this.multiply(n, c);
			r = this.multiply(r, b);
		}
		Set<PrimeFieldElement> result = new TreeSet<>();
		result.add(r);
		result.add(this.negative(r));
		return result;
	}

	@Override
	public PrimeFieldElement getRandomElement() {
		return getElement(MiscAlgorithms.randomBigInteger(new Random(), p));
	}

	@Override
	public BigInteger getNumberOfElements() {
		return this.p;
	}

	@Override
	public Iterator<PrimeFieldElement> iterator() {
		return new Iterator<PrimeFieldElement>() {
			private BigInteger i = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return this.i.compareTo(p) < 0;
			}

			@Override
			public PrimeFieldElement next() {
				BigInteger tmp = this.i;
				this.i = this.i.add(BigInteger.ONE);
				return PrimeField.this.getElement(tmp);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
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
		return new PrimeFieldElement(BigInteger.valueOf(value));
	}
	
	public PrimeFieldElement getElement(BigInteger value) {
		return new PrimeFieldElement(value);
	}

	public class PrimeFieldElement implements Element {
		private BigInteger value;

		private PrimeFieldElement(BigInteger value) {
			this.value = value.mod(p);
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(p);
		}

		public BigInteger getValue() {
			return value;
		}

		public boolean equals(Object o) {
			if (!(o instanceof PrimeFieldElement))
				return false;
			PrimeFieldElement e = (PrimeFieldElement) o;
			return this.value.equals(e.value);
		}

		public String toString() {
			return this.value.toString();
		}

		@Override
		public int compareTo(Element e) {
			PrimeFieldElement o = (PrimeFieldElement) e;
			return this.value.compareTo(o.value);
		}
	}
}
