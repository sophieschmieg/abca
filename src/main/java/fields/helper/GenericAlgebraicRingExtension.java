package fields.helper;

import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;

public class GenericAlgebraicRingExtension<T extends Element<T>>
		extends AbstractAlgebraicRingExtension<T, GenericAlgebraicExtensionElement<T>, GenericAlgebraicRingExtension<T>> {
	
	public static class GenericAlgebraicExtensionElement<T extends Element<T>>
			extends AbstractElement<GenericAlgebraicExtensionElement<T>>
			implements AlgebraicExtensionElement<T, GenericAlgebraicExtensionElement<T>> {
		private UnivariatePolynomial<T> asPolynomial;

		private GenericAlgebraicExtensionElement(UnivariatePolynomial<T> asPolynomial) {
			this.asPolynomial = asPolynomial;
		}

		@Override
		public String toString() {
			return asPolynomial.toString("α", true);
		}

		@Override
		public int compareTo(GenericAlgebraicExtensionElement<T> o) {
			return asPolynomial.compareTo(o.asPolynomial);
		}

		public UnivariatePolynomial<T> asPolynomial() {
			return asPolynomial;
		}

	}
	public GenericAlgebraicRingExtension(UnivariatePolynomial<T> minimalPolynomial, Ring<T> baseRing) {
		super(minimalPolynomial, baseRing, "x");
	}

	public GenericAlgebraicRingExtension(Ring<T> baseRing) {
		super(baseRing);
	}

	@Override
	protected GenericAlgebraicExtensionElement<T> fromSmallDegreePolynomial(UnivariatePolynomial<T> polynomial) {
		return new GenericAlgebraicExtensionElement<>(polynomial);
	}

	@Override
	public GenericAlgebraicRingExtension<T> makeExtension(UnivariatePolynomial<T> minimalPolynomial) {
		return new GenericAlgebraicRingExtension<>(minimalPolynomial, getRing());
	}

	@Override
	protected GenericAlgebraicRingExtension<T> asExtensionType() {
		return this;
	}

}
