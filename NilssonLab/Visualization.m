
(* Some useful plotting functions *)

BeginPackage["NilssonLab`Visualization`", {"NilssonLab`Utilities`"}];


SimpleBarChart::usage = "SimpleBarChart[x] plots a plain bar chart for a data vector x. SimpleBarChart[{x1, x2 ... }] plot several groups of bars.";
ErrorBar;
ErrorBarChart::usage = "ErrorBarChart[x, e] where x is a list of point and e is a corresponding list of errors draws a simple bar chart with error bars. ErrorBarChart[x, xmin, xmax] draws asymmetric error bars from xmin to xmax. ErrorBarChart[X] where X is a matrix uses the means and standard deviations of each column as errorbars.";
DataPointsPlot;
ScatterBarChart;


Begin["NilssonLab`Visualization`Private`"];


(* Drawing bar chart, in a simpler / clearer format than the built-in BarChart 
   We assume all groups are of the same size. Each group is draw to fit on a segment
   of length 1 on the x-axis. With m bars in the length 1 segment, we specify the bar width w
   as a fraction f of 1/m. We have distance d between bars within the same group,
   and distance e between bars in different groups, so that it always holds that m*w + (m-1)*d + e = 1.
   Having fixed w < 1/m, we specify the remaining parameter g = d/e to control bar spacing
   within/between groups. From  m*w + (m-1)*d + d / g = 1  we obtain
   d = g (1 - m*w)/(g(m-1) + 1) and e = 1 - m*w - (m-1)*d.
   Setting g to zero gives d = 0 si that bars within a group are touching.
   *)

SimpleBarChart[x_?VectorQ, opts:OptionsPattern[]] :=
	SimpleBarChart[Partition[x, 1], opts]

Options[SimpleBarChart] = {"FractionalBarWidth" -> 0.8, "BarSpacingRatio" -> 1/3}

SimpleBarChart[Y:{_?VectorQ..}, opts:OptionsPattern[{SimpleBarChart, Graphics}]]:=
	Block[{n, m, f, g, X, w, col},
		n = Length[Y]; (* number of groups *)
		m = Length[First[Y]]; (* no. samples per group *)
		f = OptionValue["FractionalBarWidth"];
		g = OptionValue["BarSpacingRatio"];
		X = allBarMidPoints[n, m, f, g];
		w = barWidth[m, f];
		(* distinguish by color within group i *)
		col = GrayLevel /@ (Range[0, m-1]/m);
		Graphics[
			Table[
				{col[[j]], Rectangle[{X[[i,j]] - w/2, 0}, {X[[i,j]] + w/2, Y[[i,j]]}]},
				{i, n}, {j, m}],
			FilterRules[{opts}, Options[Graphics]],
			Axes->True, AxesOrigin->{0,0}, Ticks->{None, Automatic}]]

barWidth[m_Integer, f_] := f / m

allBarMidPoints[n_Integer, m_Integer, f_, g_] :=
	Block[{w, d, e, x, col},
		w = barWidth[m, f];
		d = g * (1 - m*w) / (g*(m - 1) + 1);
		e = 1 - m*w - (m-1)*d;
		Table[
			barMidPoint[i, j, w, d, e],
			{i, n}, {j, m}]]

barMidPoint[i_Integer, j_Integer, w_, d_, e_, opts:OptionsPattern[]] := 
	(i-1) + e/2 + w/2 + (j-1)*(w + d)



(* An error bar from ymin to ymax with whiskers of width w *)

ErrorBar[x_, ymin_, ymax_, w_]:=
	{Line[{{x, ymin},{x, ymax}}],
	 Line[{{x-w/2, ymin},{x+w/2, ymin}}],
	 Line[{{x-w/2, ymax},{x+w/2, ymax}}]}

(* Draw an error bar at expected bar mid points for a list of values, with width w *)

drawErrorBars[ymin_?VectorQ, ymax_?VectorQ, f_, g_] :=
	drawErrorBars[Partition[ymin, 1], Partition[ymax, 1], f, g]

(*For grouped data. Assuming all groups are the same size*)

drawErrorBars[ymin_?MatrixQ, ymax_?MatrixQ, f_, g_] :=
	Block[{n, m, X},
		{n, m} = Dimensions[ymin];
		X = allBarMidPoints[n, m, f, g];
		w = barWidth[m, f] / 2;
		Graphics[
			Table[
				ErrorBar[X[[i,j]], ymin[[i,j]], ymax[[i,j]], w],
				{i, n}, {j, m}]]]


(* Bar chart wih error bars, for a single vector of values *)

Options[ErrorBarChart] = {"FractionalBarWidth" -> 0.8, "BarSpacingRatio" -> 1/3}

ErrorBarChart[y_?VectorQ, e_?VectorQ, opts:OptionsPattern[]]:=
	ErrorBarChart[y, y-e, y+e, opts]


(* Error bar chart from specified values and asymmetric errors, draws from xmin to xmax *)

ErrorBarChart[y_?VectorQ, ymin_?VectorQ, ymax_?VectorQ, opts:OptionsPattern[{ErrorBarChart, Graphics, Show}]]:=
	Show[
		SimpleBarChart[y, opts],
		drawErrorBars[ymin, ymax, OptionValue["FractionalBarWidth"], OptionValue["BarSpacingRatio"]],
		FilterRules[{opts}, Options[Show]]]

(* Error bar chart from a matrix, using means and standard deviations across columns *)

ErrorBarChart[data:{_?VectorQ..}, opts:OptionsPattern[]]:=
	Block[{y, err},
		y = Mean /@ data;
		err = Table[If[Length[x]==1, 0, StandardDeviation[x]], {x, data}];
		ErrorBarChart[y, err, opts]]

(* Grouped error bars, from data matrix X = {{group 1}, {group 2} ... }
   (all groups same size) and corresponding error matrix *)

ErrorBarChart[Y_?MatrixQ, Err_?MatrixQ, opts:OptionsPattern[{ErrorBarChart, Show}]]:=
	Show[
		SimpleBarChart[Y, opts],
		drawErrorBars[Y - Err, Y + Err, 
			Quiet[OptionValue["FractionalBarWidth"]],
			Quiet[OptionValue["BarSpacingRatio"]]],
		FilterRules[{opts}, Options[Show]]]

(*From a list of grouped data, taking means and standard deviations per group*)

ErrorBarChart[data_List/;Depth[data]==4, opts:OptionsPattern[]]:=
	Block[{Y, Err},
		Y = Transpose[Map[Mean, data, {2}]];
		Err = Transpose[Map[
			If[VectorQ[#]&&Length[#]>1, StandardDeviation[#], 0]&,
			data, {2}]];
		ErrorBarChart[Y, Err, opts]]


(* Plot data points ("scatter") for a series of data points, same x-axis geometry as SimpleBarChart.
  Here we have an additonal option controlling how wide the points are set relative to the bar width:
  1 = over entire bar, 0 = on top of each other along the bar midpoint *)

Options[DataPointsPlot] = {"FractionalBarWidth" -> 0.8, "BarSpacingRatio" -> 1/3, "PointSpacing" -> 0.7}

DataPointsPlot[Y:{_?VectorQ..}, opts:OptionsPattern[{DataPointsPlot, ListPlot}]] :=
	DataPointsPlot[Transpose[Partition[Y,1]], opts]


(* For a 3D data array Y with dimensions {samples x group x replicates} *)
(* Note: why not have the group dimension outermost? *)

DataPointsPlot[Y_List/;Depth[Y]==4, opts:OptionsPattern[{DataPointsPlot, ListPlot}]]:=
	Block[{n, m, f, g, X, w, xoffset},
		{m, n} = Take[Dimensions[Y], 2];
		f = OptionValue["FractionalBarWidth"];
		g = OptionValue["BarSpacingRatio"];
		X = allBarMidPoints[n, m, f, g];
		w = barWidth[m, f] * OptionValue["PointSpacing"];
		ListPlot[
			Table[
				Join @@ Table[
					xoffset = If[w > 0,
						Range[-w/2, w/2, w / (Length[Y[[j, i]]]-1)],
						Table[0, {Length[Y[[j, i]]]}]];
					Transpose[{X[[i, j]] + xoffset, Y[[j, i]]}],
					{j, m}],
				{i, n}],
			FilterRules[{opts}, Options[ListPlot]],
			Axes->True, AxesOrigin->{0,0}, Ticks->{None,Automatic}]]


(* ScatterBarChart overlays SimpleBarChart and DataPointsPlot, as alternative to ErrorBarChart *)

(*From a matrix, using means and standard deviations across columns*)

ScatterBarChart[data:{_?VectorQ..}, opts:OptionsPattern[]]:=
	Block[{y},
		y = Mean /@ data;
		Show[
			SimpleBarChart[y, opts],
			DataPointsPlot[data, opts],
			FilterRules[{opts}, Options[Show]]]]

(*From a list of grouped data, taking means over groups*)

ScatterBarChart[data_List/;Depth[data]==4, opts:OptionsPattern[]]:=
	Block[{Y},
		Y = Transpose[Map[Mean, data, {2}]];
		Show[
			SimpleBarChart[Y, opts],
			DataPointsPlot[data, opts],
			FilterRules[{opts}, Options[ListPlot]]]]


End[];
	
EndPackage[];
