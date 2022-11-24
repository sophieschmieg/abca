package fields.numberfields;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.ProductRing;
import fields.helper.ProductRing.ProductElement;
import fields.integers.Rationals;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Vector;

class NumberFieldStuffTest {

	@Test
	void test() throws IOException {
		PAdicField q2 = new PAdicField(2, 128);
		UnivariatePolynomialRing<PAdicNumber> polynomials = q2.getUnivariatePolynomialRing();
		UnivariatePolynomial<PAdicNumber> polynomial = polynomials.add(polynomials.getVarPower(2), polynomials.getInteger(7));
		Map<PAdicNumber, Integer> roots = q2.roots(polynomial);
		System.out.println(q2.roundToInteger(roots.keySet().iterator().next(), 32));
		System.out.println(q2.roundToInteger(roots.keySet().iterator().next(), 64));
			Rationals q = Rationals.q();
		NumberField nf = NumberField.getNumberField(q.getUnivariatePolynomialRing().parse("X^2 + -2"));
		int prime = 43;
		NumberFieldIdeal ideal = nf.maximalOrder().idealsOver(prime).get(0);
		ModuloNumberFieldIdeal mod = ideal.modOut();
		NFE numerator = nf.add(nf.one(), nf.alpha());
		NFE denominator = nf.subtract(nf.getInteger(2), nf.alpha());
		NFE t = mod.lift(mod.divide(mod.reduce(numerator), mod.reduce(denominator)));
		System.out.println(t);
		FreeModule<NFE> freeNf = new FreeModule<>(nf.maximalOrder(), 2);
		Vector<NFE> v1 = new Vector<>(t, nf.one());
		Vector<NFE> v2 = new Vector<>(nf.getInteger(prime), nf.zero());
		System.out.println(v1);
		System.out.println(v2);
		List<Vector<NFE>> basis = new ArrayList<>();
		basis.add(v1);
		basis.add(v2);
		FreeSubModule<NFE, Vector<NFE>> lattice = new FreeSubModule<>(freeNf, basis);
		Vector<NFE> result = new Vector<>(numerator, denominator);
		System.out.println(lattice.asVector(result));
		Reals r = Reals.r(128);
		ProductRing<Real, Reals> product = ProductRing.power(r, nf.realEmbeddings().size());
		FreeModule<ProductElement<Real>> free = new FreeModule<>(product, 2);
		List<Vector<ProductElement<Real>>> orthogonal = new ArrayList<>();
		FiniteRealVectorSpace flatSpace = new FiniteRealVectorSpace(r, 4);
		while(true) {
			orthogonal.clear();
			orthogonal.add(embed(basis.get(0), nf, product));
			orthogonal.add(embed(basis.get(1), nf, product));
				ProductElement<Real> coeff = product.divide(free.innerProduct(orthogonal.get(0), orthogonal.get(1)),
					free.innerProduct(orthogonal.get(0), orthogonal.get(0)));
			System.out.println(orthogonal);
			orthogonal.set(1, free.subtract(orthogonal.get(1), free.scalarMultiply(coeff, orthogonal.get(0))));
			System.out.println(coeff);
			NFE rounded = round(coeff, nf);
			System.out.println(rounded);
			System.out.println(orthogonal);
			basis.set(1, freeNf.subtract(basis.get(1), freeNf.scalarMultiply(rounded, basis.get(0))));
			System.out.println(basis);
			if (flatSpace.valueNorm(flatEmbed(basis.get(0), nf))
					.compareTo(flatSpace.valueNorm(flatEmbed(basis.get(1), nf))) < 0) {
				break;
			}
			Vector<NFE> tmp = basis.get(0);
			basis.set(0, basis.get(1));
			basis.set(1, tmp);
		}
	}

	private NFE round(ProductElement<Real> coeff, NumberField nf) {
		return nf.maximalOrder().getVectorSpace().closestLatticePoint(new Vector<>(coeff.values()), nf.maximalOrder());
	}

	private Vector<ProductElement<Real>> embed(Vector<NFE> t, NumberField nf, ProductRing<Real, Reals> product) {
		List<ProductElement<Real>> result = new ArrayList<>();
		for (NFE c : t.asList()) {
			List<Real> elements = new ArrayList<>();
			for (EmbeddedNumberField<Real, Reals> embedding : nf.realEmbeddings()) {
				elements.add(embedding.embedding(c));
			}
			result.add(product.getElement(elements));
		}
		return new Vector<>(result);
	}

	private Vector<Real> flatEmbed(Vector<NFE> t, NumberField nf) {
		List<Real> result = new ArrayList<>();
		for (NFE c : t.asList()) {
			result.addAll(nf.minkowskiEmbedding(c).asList());
		}
		return new Vector<>(result);
	}
}
