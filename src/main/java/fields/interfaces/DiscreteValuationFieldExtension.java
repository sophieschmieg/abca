package fields.interfaces;

import fields.local.LocalRingExtension;

public interface DiscreteValuationFieldExtension<B extends Element<B>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, FE extends DiscreteValuationFieldExtension<B, S, E, FE, R, RE, RFE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
		extends FieldExtension<B, E, FE>, DiscreteValuationField<E, RE>,
		GlobalFieldExtension<B, B, S, E, E, R, RE, RFE, FE, FE> {
	
	@Override
	DiscreteValuationField<B, S> getBaseField();
	
	@Override
	ExtensionOfDiscreteValuationField<E, RE, B, S, E, FE, R, RE, RFE> getUniqueExtension(
			UnivariatePolynomial<E> minimalPolynomial);

	@Override
	LocalRingExtension<B, S, E, FE, R, RE, RFE> ringOfIntegers();
	
	@Override
	default ExtensionOfGlobalField<E, E, RE, B, B, S, E, E, R, RE,RFE, FE, FE> getGlobalFieldExtension(
			UnivariatePolynomial<E> minimalPolynomial) {
		return getUniqueExtension(minimalPolynomial).asGlobalFieldExtension();
	}
}
