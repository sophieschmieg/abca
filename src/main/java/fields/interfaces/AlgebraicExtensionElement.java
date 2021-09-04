package fields.interfaces;

public interface AlgebraicExtensionElement<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>> extends Element<S> {
	UnivariatePolynomial<T> asPolynomial();

	String toString();
}
