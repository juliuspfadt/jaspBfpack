//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//

import QtQuick
import QtQuick.Layouts
import JASP
import JASP.Controls

Section
{
	property bool bfTy: true
	property bool iterations: false
	property bool interactions: false

	title: 	qsTr("Options")

	Group 
	{
		title: qsTr("Interaction terms")
		visible: interactions

		// VariablesForm {
		// 	preferredHeight: 300 * preferencesModel.uiScale

		// 	AvailableVariablesList
		// 	{
		// 		name: "new"
		// 		visible:false
		// 		source: "covariates"
		// 		implicitWidth: 0
		// 	}
		// 	AssignedVariablesList
		// 	{ 
		// 		name: "interactionTerms"
		// 		source: "new"
		// 		listViewType: JASP.Interaction
		// 	}
		// }

		ComponentsList 
		{
			implicitHeight: 120 * preferencesModel.uiScale
			// headerLabels: [qsTr("Include")]
			titles: [qsTr("Include")]
			name: "interactionTerms"
			rSource: "interactionSource"
			// source: "covariates"
			rowComponent: RowLayout { 
				Text {Layout.preferredWidth: 200; text: rowValue; visible: true}
				CheckBox {Layout.preferredWidth: 100; name: "includeInteractionEffect"; checked:true}
			}
		}
	}

	Group
	{
		title: qsTr("Bayes Factor")
		// Layout.rowSpan: 2

		CheckBox
		{
			name: "logScale"
			label: qsTr("On log scale")
		}

		RadioButtonGroup
		{
			visible: bfTy
			name: "bfType"
			title: qsTr("Type")
			radioButtonsOnSameRow: false
			RadioButton { value: "fractional"; label: qsTr("Fractional"); checked: true}
			RadioButton { value: "adjusted"; label: qsTr("Adjusted fractional")}
		}
	}
	Group
	{
		title: 							qsTr("Tables")

		CheckBox 
		{
			name: "specificationTable"
			text: qsTr("Specification")
		}

		CheckBox
		{
			name: 						"coefficientsTable"
			text: 						qsTr("Coefficients")

			CIField
			{
				name: 					"ciLevel"
				text: 					qsTr("Uncertainty interval")
			}
		}
	}

	Group
	{
		
		title: qsTr("Plots")
		CheckBox
		{
			name: 						"plots"
			text: 						qsTr("Manual hypotheses plots")
		}
	}

	Group
	{
		title: 							qsTr("Additional Options")

		IntegerField
		{
			visible: iterations
			name: "iterations"
			text: qsTr("No. iterations for parameter estimation")
			defaultValue: 5000
			min: 2000
			fieldWidth: 60 * preferencesModel.uiScale

		}

		DoubleField
		{
			name: 						"seed"
			text: 						qsTr("Seed")
			defaultValue: 				100
			min: 						-999999
			max: 						999999
			fieldWidth: 				60 * preferencesModel.uiScale
		}
	}
}