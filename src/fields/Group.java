package fields;

import java.math.BigInteger;
import java.util.Map;

import util.MiscAlgorithms;

public interface Group<T extends Element> extends MathSet<T> {
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
		
		Map<BigInteger, Integer> primeDecomp = MiscAlgorithms.primeDecomposition(order);
		for (BigInteger prime : primeDecomp.keySet()) {
			int power = primeDecomp.get(prime);
			while (power > 0 && power(order.divide(prime), t).equals(neutral())) {
				power--;
				order = order.divide(prime);
			}
		}
		return order;
	}
}
