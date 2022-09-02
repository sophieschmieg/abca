package fields.interfaces;

import fields.local.LocalRingExtension;

public interface DedekindRingExtension<B extends Element<B>, I extends Element<I>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, IE extends Element<IE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DFE extends DiscreteValuationFieldExtension<B, S, E, DFE, R, RE, RFE>, GFE extends GlobalFieldExtension<B, I, S, E, IE, R, RE, RFE, DFE, GFE>>
		extends DedekindRing<IE, E, RE> {
	@Override
	GFE quotientField();

	@Override
	LocalRingExtension<B, S, E, DFE, R, RE, RFE> localize(Ideal<IE> maximalIdeal);

	@Override
	DFE localizeAndQuotient(Ideal<IE> maximalIdeal);
}
