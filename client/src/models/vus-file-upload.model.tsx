import { IHPOTerm } from "../services/sample/sample.dto";

export interface IUnprocessedVus {
  locus: string;
  type: string;
  refAllele?: string | null;
  altAllele?: string | null;
}

export interface IVusGene {
  index: number;
  vus: IUnprocessedVus;
  genes: string[];
}

export interface IVusGeneSelected {
  index: number;
  gene: string;
}

export interface ISamplePhenotypeSelected {
  sampleId: string;
  phenotypesSelected: IHPOTerm[];
}
