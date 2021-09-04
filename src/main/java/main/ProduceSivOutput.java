package main;

import java.util.Map;
import java.util.TreeMap;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.FieldEmbedding;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.Monomial;

public class ProduceSivOutput {
	public static void main(String[] args) {
		new ProduceSivOutput();
	}

	public ProduceSivOutput() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		FiniteField f256 = FiniteField.getFiniteField(f2, 8);
		System.out.println(f256.minimalPolynomial());
		Map<Monomial, PFE> f2to128polymap = new TreeMap<>();
		UnivariatePolynomialRing<PFE> f2ring = f2.getUnivariatePolynomialRing();
		f2to128polymap.put(f2ring.getMonomial(new int[] { 128 }), f2.one());
		f2to128polymap.put(f2ring.getMonomial(new int[] { 127 }), f2.one());
		f2to128polymap.put(f2ring.getMonomial(new int[] { 126 }), f2.one());
		f2to128polymap.put(f2ring.getMonomial(new int[] { 121 }), f2.one());
		f2to128polymap.put(f2ring.getMonomial(new int[] { 0 }), f2.one());
		UnivariatePolynomial<PFE> f2to128poly = f2ring.getPolynomial(f2to128polymap);
		FiniteField f2to128 = FiniteField.getFiniteField(f2to128poly, f2);
		System.out.println(f2to128.minimalPolynomial());
		UnivariatePolynomialRing<FFE> f256ring = f256.getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<FFE>> fact = f256.factorization(f256ring.getEmbedding(f2to128poly, new MathMap<>() {
			@Override
			public FFE evaluate(PFE t) {
				return f256.getEmbedding(t);
			}
		}));
		for (Polynomial<FFE> factor : fact.primeFactors()) {
			System.out.println(factor);
		}
		FieldEmbedding<PFE, FFE, FiniteField> f2to128overf256 = f256.getEmbeddedExtension(f256.getUnivariatePolynomialRing().toUnivariate(fact.firstPrimeFactor())
				);
		System.out.println(f2to128overf256.minimalPolynomial());
		System.out.println(f2to128overf256.getField().minimalPolynomial());
		System.out.println(f2to128overf256.getField().alpha());
		System.out.println(f2to128overf256.getEmbeddedField().alpha());
		System.out.println(f2to128overf256.getEmbeddedAlpha());
		//System.out.println(f2to128overf256.gammaToAlphaBetaBaseChange());
		//System.out.println(f2to128overf256.alphaBetaToGammaBaseChange());
	}
}
