import { IVus } from "../../models/view-vus.model";
import {
  ISamplePhenotypeSelected,
  IVusGene,
  IVusGeneSelected,
} from "../../models/vus-file-upload.model";

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
  multipleGenesSelection?: IVusGeneSelected[];
  samplePhenotypeSelection?: ISamplePhenotypeSelected[];
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
  areRsidsRetrieved: boolean;
  isClinvarAccessed: boolean;
  vusList: IVus[];
  multipleGenes?: IVusGene[] | null;
  uniqueSampleIds?: string[] | null;
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVus[] | null;
}

export interface IGetHPOTermsRequest {
  phenotype: string;
}

export interface IHPOTerm {
  ontologyId: string;
  name: string;
}

export interface IGetHPOTermsResponse {
  isSuccess: boolean;
  hpoTerms: IHPOTerm[];
}
