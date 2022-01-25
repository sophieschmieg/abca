package fields.interfaces;

import java.math.BigInteger;
import java.util.Map;
import java.util.TreeMap;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ring.FactorizationResult;

public interface Group<T extends Element<T>> extends Monoid<T> {
	public T neutral();
	public T inverse(T t);
	public T operate(T t1, T t2);
	
	public default T power(int n, T t) {
		return power(BigInteger.valueOf(n), t);	
	}
	
	public default T power(BigInteger n, T t) {
		T result = this.neutral();
		if (n.signum() < 0) {
			n = n.negate();
			t = this.inverse(t);
		}
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.operate(result, t);
			t = this.operate(t, t);
		}
		return result;
	}
	
	public default BigInteger discreteLogarithm(T base, T power) {
		Map<T, Integer> babySteps = new TreeMap<>();
		BigInteger n = getOrder(base);
		int m = n.sqrt().add(BigInteger.ONE).intValueExact();
		T babyStep = neutral();
		for(int i = 0; i < m; i++) {
		babySteps.put(babyStep, i);
		babyStep = operate(babyStep, base);
		}
		T baseByM = power(-m, base);
		T giantStep = power;
		for(int i = 0; i < m; i++) {
			if (babySteps.containsKey(giantStep)) {
				return BigInteger.valueOf(i).multiply(BigInteger.valueOf(m)).add(BigInteger.valueOf(babySteps.get(giantStep)));
			}
			giantStep = operate(giantStep, baseByM);
		}
		throw new ArithmeticException("Could not find discrete logarithm!");
	}
	
	public default BigInteger getOrder(T t) {
		BigInteger order = this.getNumberOfElements();
		if (!power(order, t).equals(neutral())) {
			throw new ArithmeticException("Element order wrong");
		}
		
		FactorizationResult<IntE, IntE> primeDecomp = Integers.z().uniqueFactorization(new IntE(order));
		for (IntE prime : primeDecomp.primeFactors()) {
			int power = primeDecomp.multiplicity(prime);
			while (power > 0 && power(order.divide(prime.getValue()), t).equals(neutral())) {
				power--;
				order = order.divide(prime.getValue());
			}
		}
		return order;
	}
}
