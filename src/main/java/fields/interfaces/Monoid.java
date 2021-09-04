package fields.interfaces;

import java.math.BigInteger;

public interface Monoid<T extends Element<T>> extends MathSet<T> {
	public T neutral();

	public T operate(T t1, T t2);

	public default T power(BigInteger n, T t) {
		T result = this.neutral();
		if (n.signum() < 0) {
			throw new ArithmeticException("Can't take negative powers in a Monoid");
		}
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.operate(result, t);
			t = this.operate(t, t);
		}
		return result;
	}

	public default T power(int n, T t) {
		return power(BigInteger.valueOf(n), t);
	}

}
