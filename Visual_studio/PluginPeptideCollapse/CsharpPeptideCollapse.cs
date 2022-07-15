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
        public override string Name => "Peptide collapse v1.4.4";
        public override string Description => "Collapse peptides with shared sequence but different modifications into consensus sequence. Further allows stoichiometry calculation.";
        
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
            

            var ConditionColumn1Param = new SingleChoiceParam("Condition column, eg R.FileName")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("r.filename"))),
                Help = "Chose experimental condition column, such as R.FileName, which will be used to create wide-format intensity columns. Has to be a text column."
            };

            var PerlParam = new SingleChoiceWithSubParams("Condition grouping")
            {
                Values = new[] { "Group by condition column", "Skip" },
                Value = 0,
                Help = "Chose if input should be grouped from long-format (e.g. Spectronaut report or MaxQuant LFQ evidence) or mixed-format " +
                "(e.g. MaxQuant TMT evidence with conditions) into wide-format required for this plugin. Chose skip if data is already in wide-format " +
                "(e.g. MaxQuant TMT evidence without conditions).",
                SubParams = new Parameters[] { // Define the different sub-parameter groups
				   new Parameters(ConditionColumn1Param), // condition grouping
				   new Parameters() // Skip
               }
            };



            var Cutoff1aParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff1bParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff2aParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff2bParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff2cParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff3aParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };
            var Cutoff3bParam = new DoubleParam("Localization cutoff", 0.75)
            {
                Help = "Localization cutoff to be used to filter PTM localization. Can be between 0 and 1. Setting to 0 will not filter anything."
            };

            var PTMPosCol1Param = new SingleChoiceParam("PTM position column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PrecursorId (SN)" },
                Help = "Chose which column PTM positions should be extracted from. Requires respective columns."
            };
            var PTMPosCol2Param = new SingleChoiceParam("PTM position column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PrecursorId (SN)", "Modified sequence (MQ)" },
                Help = "Chose which column PTM positions should be extracted from. Requires respective columns."
            };
            var PTMPosCol3Param = new SingleChoiceParam("PTM position column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PrecursorId (SN)" },
                Help = "Chose which column PTM positions should be extracted from. Requires respective columns."
            };

            var PTMPosCol2cParam = new SingleChoiceParam("Target PTM Probabilities column")
            {
                Values = data.StringColumnNames,
                Value = Math.Max(0, data.StringColumnNames.FindIndex(col => col.ToLower().Equals("phospho (sty) probabilities"))),
                Help = "Chose which column PTM target probabilities should be extracted from. Requires respective columns."
            };

            var ProbType1Param = new SingleChoiceWithSubParams("Probability column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PTMAssayProbability (SN)", "No probability column" },
                Help = "Chose which column probabilities and PTM positions should be extracted from. Requires resepctive columns. " +
                "EG.PTMAssayProbability will extract PTM positions from EG.PrecursorId.",
                SubParams = new Parameters[] { // Define the different sub-parameter groups
				   new Parameters(Cutoff1aParam), // EG.PTMLocalizationProbabilities
				   new Parameters(Cutoff1bParam), // EG.PTMAssayProbability
                   new Parameters(PTMPosCol1Param) // Ignored
			   }
            };
            var ProbType2Param = new SingleChoiceWithSubParams("Probability column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PTMAssayProbability (SN)",
                    "No probability column", "Modified sequence (MQ)" },
                Help = "Chose which column probabilities and PTM positions should be extracted from. Requires resepctive columns. " +
                "EG.PTMAssayProbability will extract PTM positions from EG.PrecursorId. Modified sequence will extract probabilities from provided PTM Probabilities columns.",
                SubParams = new Parameters[] { // Define the different sub-parameter groups
				   new Parameters(Cutoff2aParam), // EG.PTMLocalizationProbabilities
				   new Parameters(Cutoff2bParam), // EG.PTMAssayProbability
                   new Parameters(PTMPosCol2Param), // Ignored
                   new Parameters(Cutoff2cParam, PTMPosCol2cParam) // Modified sequence
			   }
            };
            var ProbType3Param = new SingleChoiceWithSubParams("Probability column type")
            {
                Values = new[] { "EG.PTMLocalizationProbabilities (SN)", "EG.PTMAssayProbability (SN)", "No probability column" },
                Help = "Chose which column probabilities and PTM positions should be extracted from. Requires resepctive columns. " +
                "EG.PTMAssayProbability will extract PTM positions from EG.PrecursorId.",
                SubParams = new Parameters[] { // Define the different sub-parameter groups
				   new Parameters(Cutoff3aParam), // EG.PTMLocalizationProbabilities
				   new Parameters(Cutoff3bParam), // EG.PTMAssayProbability
                   new Parameters(PTMPosCol3Param) // Ignored
			   }
            };

            var GenesParam = new SingleChoiceParam("Genes or protein groups")
            {
                Values = new[] { "PG.Genes", "PG.ProteinGroups" },
                Help = "Chose if sites should be collapsed on gene or protein level."
            };
            var FastaFileParam = new FileParam("FASTA file (optional)")
            {
                Help = "Define file location of FASTA file for PTM motif sequence annotation. Leaving blank will skip annotation."
            };
            var FastaStringParam = new StringParam("FASTA identifier rule", ".*GN=([^ ]*) .*")
            {
                Help = "If FASTA file is provided, define identifier parsing rule. This is crucial to allow gene or protein name matching as defined above. Set e.g." +
                ".*GN=([^ ]*) .* for gene readout, or .*\\|(.*)\\|.* for protein readout"
            };
            var StoichParam = new SingleChoiceParam("Stoichiometry calculation")
            {
                Values = new[] { "Calculate stoichiometries", "Skip" },
                Help = "Chose if target PTM stoichiometries should be calculated."
            };


            var CollapseParam = new SingleChoiceWithSubParams("Collapse level")
            {
                Values = new[] { "Target PTM site-level, e.g. ABC1_S15_M1", "Target PTM peptide-level, e.g. FS[ph]EAMST[ph]R (stoichiometry possible)",
                    "ModSpec peptide-level, e.g. FSEAMSTR_2[ph];1[ox]" },
                Value = 0,
                Help = "Chose if precursors should be collapsed on peptide- or site-level. Site-level requires PEP.PeptidePosition column.",
                SubParams = new Parameters[] { // Define the different sub-parameter groups
                   new Parameters(ProbType1Param, GenesParam, FastaFileParam, FastaStringParam), // Localized Site
                   new Parameters(ProbType2Param, StoichParam), // Localized Peptide
				   new Parameters(ProbType3Param) // ModSpec Peptide
			   }
            };



            var PTMTypesParam = new StringParam("Variable PTMs, target PTM first", "[Phospho (STY)];[Deamidation (NQ)];[Oxidation (M)];[Carbamidomethyl (C)]")
            {
                Help = "List PTMs as listed in EG.PrecursorId, separated by semicolon. If target site- or target peptide-level collapse is performed, " +
                "target PTM needs to be listed first."
            };

            

            var AggParam = new SingleChoiceParam("Aggregation type")
            {
                Values = new[] { "Linear modeling based", "Summing" },
                Help = "Chose if peptide intensities should be aggregated based on a linear model or simply summed."
            };


            var NumberOfCoresParam = new IntParam("CPUcores", 8)
            {
                Help = "Number of CPU threads to be created by the plugin where possible. Set to 1 for lowest CPU usage. Do not set higher than actual number of available CPU threads."
            };


            return new Parameter[]
            {
                PerlParam,
                CollapseParam,
                PTMTypesParam,
                AggParam,
                NumberOfCoresParam,
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
