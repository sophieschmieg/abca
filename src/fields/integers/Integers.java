package fields.integers;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.finitefields.PrimeField.PrimeFieldElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Ideal;

public class Integers extends AbstractRing<IntE> {
	private static Integers z = new Integers();

	private Integers() {
	}

	public static Integers z() {
		return z;
	}

	public static class IntE implements Element<IntE> {
		private BigInteger value;

		public IntE(int value) {
			this.value = BigInteger.valueOf(value);
		}

		public IntE(BigInteger value) {
			this.value = value;
		}

		@Override
		public int compareTo(IntE other) {
			return this.value.compareTo(((IntE) other).value);
		}

		public boolean equals(Object o) {
			return this.value.equals(((IntE) o).value);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}
	}

	@Override
	public IntE getRandomElement() {
		return new IntE((int) Math.round(new Random().nextGaussian() * 10.0));
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<IntE> iterator() {
		throw new InfinityException();
	}

	@Override
	public IntE zero() {
		return new IntE(0);
	}

	@Override
	public IntE one() {
		return new IntE(1);
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public IntE add(IntE t1, IntE t2) {
		return new IntE(t1.value.add(t2.value));
	}

	@Override
	public IntE negative(IntE t) {
		return new IntE(t.value.negate());
	}

	@Override
	public IntE multiply(IntE t1, IntE t2) {
		return new IntE(t1.value.multiply(t2.value));
	}

	@Override
	public IntE getInteger(BigInteger t) {
		return new IntE(t);
	}

	@Override
	public boolean isUnit(IntE t) {
		return t.value.equals(BigInteger.ONE) || t.value.equals(BigInteger.ONE.negate());
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return BigInteger.TWO;
	}

	@Override
	public Iterable<IntE> getUnits() {
		List<IntE> list = new ArrayList<>();
		list.add(new IntE(1));
		list.add(new IntE(-1));
		return list;
	}

	@Override
	public IntE inverse(IntE t) {
		return t;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(IntE t) {
		return t.value.equals(BigInteger.ZERO);
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isDivisible(IntE dividend, IntE divisor) {
		return dividend.value.mod(divisor.value).equals(BigInteger.ZERO);
	}

	@Override
	public BigInteger euclidMeasure(IntE t) {
		return t.value.abs();
	}

	@Override
	public List<IntE> quotientAndRemainder(IntE dividend, IntE divisor) {
		BigInteger[] result = dividend.value.divideAndRemainder(divisor.value);
		List<IntE> list = new ArrayList<>();
		list.add(new IntE(result[0]));
		list.add(new IntE(result[1]));
		return list;
	}
	
	@Override
	public IntE projectToUnit(IntE t) {
		int signum = t.value.signum();
		if (signum == 0) {
			return one();
		}
		return getInteger(signum);
	}

	@Override
	public int krullDimension() {
		return 1;
	}
	
	@Override
	public String toString() {
		return "Z";
	}
	
	@Override
	public Ideal<IntE> getIdeal(List<IntE> generators) {
		if (generators.size() == 0) {
			return new IntegerIdeal(0);
		}
		IntE generator = generators.get(0);
		for (int i = 1; i < generators.size(); i++) {
			generator = gcd(generator, generators.get(i));
		}
		return new IntegerIdeal(generator);
	}

	@Override
	public Ideal<IntE> intersect(Ideal<IntE> t1, Ideal<IntE> t2) {
		return multiply(t1, t2);
	}

	public IntE lift(ModularIntegerRingElement t) {
		return new IntE(t.getValue());
	}

	public IntE lift(PrimeFieldElement t) {
		return new IntE(t.getValue());
	}

	public static class IntegerIdeal extends AbstractIdeal<IntE> {
		private IntE m;

		public IntegerIdeal(IntE m) {
			super(z);
			this.m = m;
		}

		public IntegerIdeal(int m) {
			this(BigInteger.valueOf(m));
		}

		public IntegerIdeal(BigInteger m) {
			this(new IntE(m));
		}

		@Override
		public boolean isFinite() {
			return false;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new InfinityException();
		}

		@Override
		public List<IntE> generators() {
			return Collections.singletonList(m);
		}

		@Override
		public List<IntE> generate(IntE t) {
			return Collections.singletonList(z.getInteger(t.value.divide(m.value)));
		}

		@Override
		public IntE residue(IntE t) {
			return z.getInteger(t.value.mod(m.value));
		}
		
		@Override
		public boolean isPrime() {
			return m.value.isProbablePrime(100);
		}

		@Override
		public boolean isMaximal() {
			return isPrime() && !m.equals(z.zero());
		}

		@Override
		public boolean contains(IntE t) {
			return t.value.mod(m.value).equals(BigInteger.ZERO);
		}
	}
}
