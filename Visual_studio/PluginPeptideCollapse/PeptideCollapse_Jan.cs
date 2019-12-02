using System;
using System.IO;
using System.Text;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Matrix;
using PluginInterop;
using Path = System.IO.Path;
using Utils = PluginInterop.R.Utils;

namespace PluginPeptideCollapse
{
	public class CsharpPeptideCollapse : PluginInterop.R.MatrixProcessing
	{
		public override string Name => "Peptide collapse v1.1";
		public override string Description => "Collapse peptides with shared sequence but different modifications into consensus sequence";

		protected override bool TryGetCodeFile(Parameters param, out string codeFile)
		{
			byte[] code = (byte[])global::PluginPeptideCollapse.Properties.Resources.ResourceManager.GetObject("RcodePeptideCollapse");
			codeFile = Path.GetTempFileName();
			File.WriteAllText(codeFile, Encoding.UTF8.GetString(code));
			return true;
		}

		protected override Parameter[] SpecificParameters(IMatrixData data, ref string errString)
		{
			if (data.StringColumnCount == 0)
			{
				errString = "Please add at least one string column";
				return null;
			}
            if (data.CategoryColumnCount == 0)
            {
                errString = "Please add at least one categorical column";
                return null;
            }
            var NumberOfCoresParam = new StringParam("CPUcores", "8")
            {
                Help = "Number of CPU cores to be used by the plugin in certain steps, set to 1 for lowest CPU usage."
            };
            var ConditionColumnParam = new SingleChoiceParam("Condition column, eg R.FileName")
            {
                Values = data.CategoryColumnNames,
                Value = Math.Max(0, data.CategoryColumnNames.FindIndex(col => col.ToLower().Equals("r.filename"))),
                Help = "Chose experimental condition column, such as R.FileName. Has to be a categorical column."
            };
            var ConditionOldParam = new StringParam("Rename condition old", "20171125_QE7_nLC14_DBJ_SA_DIAphos_RPE1_pilot2_")
            {
                Help = "Write parsing rule of condition column to be changed from. Leave blank if no change wanted."
            };
            var ConditionNewParam = new StringParam("Rename condition new", "DIA_")
            {
                Help = "Write parsing rule of condition column to be changed into. Leave blank if no change wanted."
            };
            var PerlParam = new SingleChoiceParam("Parsing type")
            {
                Values = new[] { "fixed", "perl" },
                Help = "Chose if perl compatibility should be enabled for parsing."
            };
	    var ProbParam = new SingleChoiceWithSubParams("Prob param") // Hierarchical parameter with sub parameters
	    {
		    Values = new [] {"Local", "Global"}, // Name the different sub-parameter groups
			   SubParams = new Parameters[] { // Define the different sub-parameter groups
				   new Parameters(ConditionColumnParam, ConditionOldParam, ConditionNewParam, PerlParam), // Local
				   new Parameters(PTMProbParam) // Global
			   }
	    };
            var PTMProbParam = new SingleChoiceParam("EG.PTMLocalizationProbabilities")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("eg.ptmlocalizationprobabilities"))),
                Help = "Chose localization probability column, such as EG.PTMLocalizationProbabilities. Has to be a string column."
            };
            var PTMSeqParam = new SingleChoiceParam("EG.PrecursorId")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("eg.precursorid"))),
                Help = "Chose sequence column, such as EG.PrecursorId. Has to be a string column."
            };
            var GenesParam = new SingleChoiceParam("PG.Genes or PG.ProteinGroups")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("pg.genes"))),
                Help = "Chose gene or protein group column. Has to be a string column."
            };
            var PeptidePosParam = new SingleChoiceParam("PEP.PeptidePosition")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("pep.peptideposition"))),
                Help = "Chose peptide position column, such as PEP.PeptidePosition. Has to be a string column."
            };
            var AggParam = new SingleChoiceParam("Aggregation type")
            {
                Values = new[] { "Linear modeling based", "Summing" },
                Help = "Chose if peptide intensities should be aggregated based on a linear model or simply summed."
            };
            var CutoffParam = new StringParam("Localization cutoff", "0")
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            return new Parameter[]
            {
                NumberOfCoresParam,
		ProbParam, // Use ProbParam instead of other parameters.
                PTMSeqParam,
                GenesParam,
                PeptidePosParam,
                AggParam,
                CutoffParam,
            };
        }

		protected override string GetCommandLineArguments(Parameters param)
		{
			var tempFile = Path.GetTempFileName();
			param.ToFile(tempFile);
			return tempFile;
		}

		protected override string CodeFilter => "R script, *.R | *.R";
		protected override FileParam ExecutableParam()
		{
			return Utils.CreateCheckedFileParam(InterpreterLabel, InterpreterFilter, TryFindExecutable);
		}

		protected override bool TryFindExecutable(out string path)
		{
			return Utils.TryFindRExecutable(out path);
		}
	}
}
