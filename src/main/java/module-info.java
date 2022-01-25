module abca {
	opens fields.finitefields;
	opens fields.floatingpoint;
	opens fields.integers;
	opens fields.interfaces;
	opens fields.local;
	opens fields.numberfields;
	opens fields.polynomials;
	opens fields.vectors;
	opens util;
	opens varieties.curves;
	requires java.desktop;
}