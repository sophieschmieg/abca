package fields.interfaces;

import java.math.BigInteger;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ring.FactorizationResult;

public interface Group<T extends Element<T>> extends Monoid<T> {
	public T neutral();
	public T inverse(T t);
	public T operate(T t1, T t2);
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
