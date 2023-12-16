import QtQuick 		2.12
import JASP.Module 	1.0

Description
{
	name			: "jaspBfpack"
	title			: qsTr("BFpack")
	description		: qsTr("A module for computing Bayes factors for equality, inequality, and order constrained hypotheses.")
	icon			: "bain-module"
	requiresData	: true
	version			: "0.18.1"
	author			: "Mulder, J., Williams, D. R., Gu, X., Tomarken, A., Böing-Messing, F., Olsson-Collentine, A., Meijerink, M., Menke, J., Fox, J.-P., Hoijtink, H., Rosseel, Y., Wagenmakers, E.J., and van Lissa, C."
	maintainer		: "Julius M. Pfadt <julius.pfadt@gmail.com>"
	website			: "https://bfpack.info"
	license			: "GPL (>= 3)"
	hasWrappers		: false

	GroupTitle
	{
		title: 		"T-Tests"
		icon: 		"analysis-bain-ttest.svg"
	}
	Analysis
	{
		menu:			"Independent Samples T-Test"
		title:		"BFpack Independent Samples T-Test"
		func:			"bfpackTTestIndependentSamples"
	}
	Analysis
	{
		menu: 		"Paired Samples T-Test"
		title:		"BFpack Paired Samples T-Test"
		func:			"bfpackTTestPairedSamples"
	}
	Analysis
	{
		menu:  		"One Sample T-Test"
		title:		"BFpack One Sample T-Test"
		func:			"bfpackTTestOneSample"
	}

	GroupTitle
	{
		title: 		"ANOVA"
		icon: 		"analysis-bain-anova.svg"
	}
	Analysis
	{
		menu:   	"ANOVA"
		title:		"BFpack ANOVA"
		func:			"bfpackAnova"
	}
	Analysis
	{
		menu:   	"ANCOVA"
		title:		"BFpack ANCOVA"
		func:			"bfpackAncova"
	}

	GroupTitle
	{
		title: 		"Regression"
		icon: 		"analysis-bain-regression.svg"
	}
	Analysis
	{
		menu:   	"Linear Regression"
		title:		"Bfpack Regression"
		func:			"bfpackRegressionLinear"
	}

	GroupTitle
	{
		title: 		"Correlation"
		// icon: 		"analysis-bain-regression.svg"
	}
	Analysis
	{
		menu: 		"Correlation"
		title:		"Bfpack Correlation"
		func:			"bfpackCorrelation"
	}

	GroupTitle
	{
		title: 		"Variances"
		// icon: 		"analysis-bain-regression.svg"
	}
	Analysis
	{
		menu: 		"Variances"
		title:		"Bfpack Variances"
		func:			"bfpackVariances"
	}
}
