package fields.interfaces;

public interface GlobalFieldExtension<B extends Element<B>, I extends Element<I>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, IE extends Element<IE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DFE extends DiscreteValuationFieldExtension<B, S, E, DFE, R, RE, RFE>, FE extends GlobalFieldExtension<B, I, S, E, IE, R, RE,RFE, DFE, FE>>
		extends FieldExtension<B, E, FE>, GlobalField<E, IE, RE> {
	@Override
	DedekindRingExtension<B, I, S, E, IE,R, RE, RFE, DFE, FE> ringOfIntegers();
	
	@Override
	GlobalField<B, I, S> getBaseField();

	@Override
	ExtensionOfGlobalField<E, IE, RE, B, I, S, E, IE, R, RE, RFE, DFE, FE> getGlobalFieldExtension(
			UnivariatePolynomial<E> minimalPolynomial);
}
