package fields.numberfields;

import fields.finitefields.FiniteField.FFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.FieldExtension;
import fields.interfaces.UnivariatePolynomial;
import fields.local.PAdicField.PAdicNumber;
import fields.numberfields.CompletedNumberField.Ext;

public class CompletedNumberField extends AbstractFieldExtension<PAdicNumber, Ext, CompletedNumberField>
		implements FieldExtension<PAdicNumber, Ext, CompletedNumberField>, DiscreteValuationField<Ext, FFE> {
	public static class Ext extends AbstractElement<Ext> implements AlgebraicExtensionElement<PAdicNumber, Ext> {
		private UnivariatePolynomial<PAdicNumber> asPolynomial;

		private Ext(UnivariatePolynomial<PAdicNumber> asPolynomial) {
			this.asPolynomial = asPolynomial;
		}

		@Override
		public int compareTo(Ext o) {
			return asPolynomial.compareTo(o.asPolynomial);
		}

		@Override
		public UnivariatePolynomial<PAdicNumber> asPolynomial() {
			return asPolynomial;
		}

	}
}
