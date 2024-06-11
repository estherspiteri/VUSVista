import {
  IAddPhenotypeRequest,
  IAddPhenotypeResponse,
  IDeleteSampleRequest,
  IDeleteSampleResponse,
  IGetHPOTermsRequest,
  IGetHPOTermsResponse,
  IGetSampleRequest,
  IGetSampleResponse,
  ILoadAllSamplesResponse,
  IRemovePhenotypeRequest,
  IRemovePhenotypeResponse,
  IUpdateHgvsRequest,
  IUpdateHgvsResponse,
} from "./sample.dto";

export class SampleService {
  async loadAllSamples(): Promise<ILoadAllSamplesResponse> {
    const result: ILoadAllSamplesResponse = await fetch(`/sample/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getSample(input: IGetSampleRequest): Promise<IGetSampleResponse> {
    const result: IGetSampleResponse = await fetch(
      `/sample/view/${input.sampleId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getHPOTerms(input: IGetHPOTermsRequest): Promise<IGetHPOTermsResponse> {
    const result: IGetHPOTermsResponse = await fetch(
      `/sample/phenotype/${input.phenotype}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addPhenotype(
    input: IAddPhenotypeRequest
  ): Promise<IAddPhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await fetch(`/sample/add-phenotype`, {
      method: "POST",
      body: data,
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removePhenotype(
    input: IRemovePhenotypeRequest
  ): Promise<IRemovePhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await fetch(
      `/sample/remove-phenotype`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async deleteSample(
    input: IDeleteSampleRequest
  ): Promise<IDeleteSampleResponse> {
    const result: IGetHPOTermsResponse = await fetch(
      `/sample/delete/${input.sampleId}`,
      {
        method: "DELETE",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async updateHgvs(input: IUpdateHgvsRequest): Promise<IUpdateHgvsResponse> {
    const result: IUpdateHgvsResponse = await fetch(
      `/sample/edit-hgvs/${input.sampleId}/${input.variantId}/${input.updatedHgvs}`,
      {
        method: "POST",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const samplesService = new SampleService();
