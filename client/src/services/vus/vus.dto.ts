import { IVus } from "../../models/view-vus.model";
import {
  ISamplePhenotypeSelected,
  IVusGene,
  IVusGeneSelected,
} from "../../models/vus-file-upload.model";
import { IVUSSummary } from "../../models/vus-summary.model";
import { IVusUpload } from "../../models/vus-upload.model";

export interface IUploadVusRequest {
  vus: IVusUpload;
}

export interface IUploadVusResponse {
  isSuccess: boolean;
  vusList: IVus[];
}

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
  multipleGenesSelection?: IVusGeneSelected[];
  samplePhenotypeSelection?: ISamplePhenotypeSelected[];
}

export interface INoHPOTermPhenotype {
  phenotype: string;
  samples: string[];
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
  areRsidsRetrieved: boolean;
  isClinvarAccessed: boolean;
  vusList: IVus[];
  multipleGenes?: IVusGene[] | null;
  noHpoTermPhenotypes: INoHPOTermPhenotype[];
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVUSSummary[] | null;
}

export interface IGetVusRequest {
  vusId: number;
}

export interface IGetVusResponse {
  isSuccess: boolean;
  vus: IVus;
}

export interface IVerifyGeneRequest {
  geneName: string;
}

export interface IVerifyGeneResponse {
  isSuccess: boolean;
  geneId?: number | null;
}
